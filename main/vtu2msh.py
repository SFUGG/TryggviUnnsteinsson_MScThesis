import os
import gmsh
import math
import numpy as np
import pyvista as pv
import argparse

# Initiate the argument parser
parser = argparse.ArgumentParser()
# Add arguments
parser.add_argument('-g', '--geometry', action='store', type=str, required=True, choices=['cave', 'chimney'],
                    help='Select glaciovolcanic void geometry: "cave" or "chimney".')
parser.add_argument('-H', '--thickness', action='store', type=float, default=100,
                    help='Set glacier thickness.')
parser.add_argument('-I', '--incline', action='store', type=float, required=True,
                    help='The bed incline.')
parser.add_argument('-f', '--file', action='store', type=str, required=True,
                    help='Specify the results file that is to be remeshed.')
parser.add_argument('-m', '--minRes', action='store', type=float, default=2,
                    help='Set the minimum mesh resolution around the glaciovolcanic void.')
parser.add_argument('-M', '--maxRes', action='store', type=float, default=10,
                    help='Set the maximum mesh resolution around the boundaries.')

# Read arguments from the command line
args = parser.parse_args()

# Assign some global variables from the command line.
geometry = args.geometry        # The void geometry.
H        = args.thickness       # The initial ice thickness.
alpha    = args.incline         # The bed incline.
lcv      = args.minRes          # The mesh resolution at the void.
lc       = args.maxRes          # The mesh resolution at the boundary.
tol      = 1e-1     			# The vertical tolerance between bed and ice.

# The domain limits and centre point C.
x_min, x_max, y_min, y_max = [0, 2*H, 0, 2*H]
C = np.array([(x_max - x_min)/2, (y_max - y_min)/2, np.tan(alpha * np.pi / 180)*(x_max - x_min)/2])

# The current working directory (cwd).
cwd = os.getcwd()
# Check if the temporary folder exists and create it if not.
temp_folder = os.path.join(cwd, 'temp')
if not os.path.exists(temp_folder):
    os.mkdir(temp_folder)


def vtu2stl(vtu_file):
    """
    Extract the 6 different surfaces from the .vtu output file from Elmer/Ice, and preserve their geometrical ID's.
    Input: The .vtu file location.
    Output: The centre line of the void and a pseudo radial geometry.
    """

    # Use pyvista to read in the mesh.
    vtu = pv.read(vtu_file)

    # Extract the surfaces.
    for i in [101, 106]:
        # Find the indices of the cells that belong to surface 100+i
        ids = np.arange(vtu.n_cells)[vtu['GeometryIds'] == i]
        surf = vtu.extract_cells(ids)

        # Save the individual surface as a .stl file.
        name = os.path.join(temp_folder, 'surf' + str(i) + '.stl')
        # Save the surface again.
        pv.save_meshio(name, surf)

    return


def remesh_surface_gmsh(filename, mesh_resolution):
    """
    Remesh a surface.
    """
    # Initialize gmsh module
    gmsh.initialize()
    gmsh.option.setNumber("Mesh.Algorithm", 5)
    gmsh.option.setNumber('General.Verbosity', 1)

    # Set the mesh grid size.
    gmsh.option.setNumber("Mesh.MeshSizeMin", mesh_resolution)
    gmsh.option.setNumber("Mesh.MeshSizeMax", mesh_resolution)

    # load an STL surface
    gmsh.merge(filename)

    # Classify the surface mesh and create geometries for the surface and the curves.
    gmsh.model.mesh.classifySurfaces(math.pi, curveAngle=math.pi / 3)
    gmsh.model.mesh.createGeometry()
    gmsh.model.geo.synchronize()

    # To generate the 2D mesh
    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.optimize('', True)

    # Overwrite the .stl file with the new mesh.
    gmsh.write(filename)
    gmsh.finalize()


def align_base_points(void_points):
    """
    Fix irregularities in the void mesh at the bed boundary that result from the extraction process.
    :param void_points: The points belonging to the void, in a format of a Nx3 array.
    :return: Returns the points once the have been fixed.
    """
    # Find the points at the base.
    ind_base = np.argwhere(void_points[:, 2] - np.tan(alpha * np.pi / 180) * (x_max - void_points[:, 0]) < tol).T
    base_points = void_points[ind_base[0], :]

    # Project the base points down to the xy-plane.
    base_points_2D = np.vstack((np.sqrt(base_points[:, 2] ** 2 + (x_max - base_points[:, 0]) ** 2),
                                base_points[:, 1])).T
    base_centre_2D = np.array([base_points_2D[:, 0].mean(), C[1]])

    # Compute the distance between each point and the centre.
    base_distance = np.array([np.linalg.norm(base_point - base_centre_2D) for base_point in base_points_2D])

    # Find the indices based on the angle of the points.
    base_angle = np.angle(base_points_2D[:, 0]-base_centre_2D[0] + (base_points_2D[:, 1]-base_centre_2D[1]) * 1j)
    ind_angle = np.argsort(base_angle)
    # Save the correctly arranged array of the base points.
    base_points_new = base_points[ind_angle]

    # Fix points that are too close to the centre compared to their neighbours.
    def fix_distance(r_point, angle_point, r_pre, r_aft):
        """
        If a point protrudes from the boundary it is set as the average of its neighbours.
        """
        # Check if a point is too far in compared to its neighbours.
        if (r_pre - r_point > 0) and (r_aft - r_point > 0):
            r_point = (r_pre + r_aft) / 2
        # Check if a point is too far out compared to its neighbours.
        elif (r_point - r_pre > 0) and (r_point - r_aft > 0):
            r_point = (r_pre + r_aft) / 2

        # Adjust the 2D point accordingly.
        new_point_2D = base_centre_2D + r_point * np.array([np.cos(angle_point), np.sin(angle_point)])

        # Project the 2D point back to the slope.
        new_point = np.array([(x_max - np.cos(alpha * np.pi / 180) * new_point_2D[0]), new_point_2D[1],
                              np.sin(alpha * np.pi / 180) * new_point_2D[0]])

        return new_point

    # Fix points that have make the boundary too sharp at that point.
    def fix_curvature(r_point, angle_point, r_pre, r_aft, curve_pre, curve_aft):
        """
        If a point marks a sharp discontinuity in the boundary it is set as the average of its neighbours.
        """
        # Check if the curvature makes a drastic change (more than 15°).
        if abs(curve_pre - curve_aft) > 15 * np.pi/180:
            r_point = (r_pre + r_aft) / 2

        # Adjust the 2D point accordingly.
        new_point_2D = base_centre_2D + r_point * np.array([np.cos(angle_point), np.sin(angle_point)])

        # Project the 2D point back to the slope.
        new_point = np.array([(x_max - np.cos(alpha * np.pi / 180) * new_point_2D[0]), new_point_2D[1],
                              np.sin(alpha * np.pi / 180) * new_point_2D[0]])

        return new_point

    # Apply the fixes to all the points at the base.
    for i in range(len(base_points_new)):
        # Variables that we need for the fixes are previous, current and next radial distance, current point angle,
        #   the angle between the previous and current and the current and next point.
        var_r_point = base_distance[ind_angle][i]
        var_angle_point = base_angle[ind_angle][i]
        if i == 0:
            var_r_pre = base_distance[ind_angle][-1]
            var_r_aft = base_distance[ind_angle][i + 1]
            var_curv_pre = np.angle(base_points_new[i, 0]-base_points_new[-1, 0] +
                                    1j*(base_points_new[i, 1]-base_points_new[-1, 1]))
            var_curv_aft = np.angle(base_points_new[i + 1, 0] - base_points_new[i, 0] +
                                    1j * (base_points_new[i + 1, 1] - base_points_new[i, 1]))
        elif i == len(base_points_new) - 1:
            var_r_pre = base_distance[ind_angle][i - 1]
            var_r_aft = base_distance[ind_angle][0]
            var_curv_pre = np.angle(base_points_new[i, 0]-base_points_new[i - 1, 0] +
                                    1j*(base_points_new[i, 1]-base_points_new[i - 1, 1]))
            var_curv_aft = np.angle(base_points_new[0, 0] - base_points_new[i, 0] +
                                    1j * (base_points_new[0, 1] - base_points_new[i, 1]))
        else:
            var_r_pre = base_distance[ind_angle][i - 1]
            var_r_aft = base_distance[ind_angle][i + 1]
            var_curv_pre = np.angle(base_points_new[i, 0]-base_points_new[i - 1, 0] +
                                    1j*(base_points_new[i, 1]-base_points_new[i - 1, 1]))
            var_curv_aft = np.angle(base_points_new[i + 1, 0] - base_points_new[i, 0] +
                                    1j * (base_points_new[i + 1, 1] - base_points_new[i, 1]))

        # Compute and apply the fixes.
        base_points_new[i, :] = fix_distance(r_point=var_r_point,
                                             angle_point=var_angle_point,
                                             r_pre=var_r_pre,
                                             r_aft=var_r_aft)
        base_points_new[i, :] = fix_curvature(r_point=var_r_point,
                                              angle_point=var_angle_point,
                                              r_pre=var_r_pre,
                                              r_aft=var_r_aft,
                                              curve_pre=var_curv_pre,
                                              curve_aft=var_curv_aft)

        # Re-compute the distance between each point and the centre.
        base_distance = np.array([np.linalg.norm(base_point - base_centre_2D) for base_point in base_points_2D])

        # Find the indices based on the angle of the points.
        base_angle = np.angle(base_points_2D[:, 0] - base_centre_2D[0] +
                              (base_points_2D[:, 1] - base_centre_2D[1]) * 1j)

    # Transfer the updated points to the actual base points.
    base_points[ind_angle] = base_points_new

    # Update the total array void points.
    void_points[ind_base, :] = base_points

    return void_points


def extract_and_remesh_void():
    """
    Extract the void from the bottom surface (101).
    """
    # The bottom mesh.
    stl_file = os.path.join(temp_folder, 'surf101.stl')
    # Use pyvista to read in the mesh.
    stl = pv.read(filename=stl_file)

    # Make sure that points below the bed or within the prescribed tolerance line up with the bed.
    ind_lower = stl.points[:, 2] - np.tan(alpha * np.pi / 180) * (x_max - stl.points[:, 0]) < 0
    stl.points[:, 2][ind_lower] = np.tan(alpha * np.pi / 180) * (x_max - stl.points[:, 0][ind_lower])
    ind_bed = np.abs(stl.points[:, 2] - np.tan(alpha * np.pi / 180) * (x_max - stl.points[:, 0])) <= tol
    stl.points[:, 2][ind_bed] = np.tan(alpha * np.pi / 180) * (x_max - stl.points[:, 0][ind_bed])

    # Find the indices of the cells that belong to the void, based on the prescribed tolerance.
    i_tf = [(abs(stl.cell_points(i)[:, 2] -
                 np.tan(alpha * np.pi / 180) * (x_max - stl.cell_points(i)[:, 0])) > tol).any()
            for i in range(stl.n_cells)]
    ind = np.arange(stl.n_cells)[i_tf]
    # Extract the void.
    void = stl.extract_cells(ind=ind)

    # Fix the void points at the base and do so until no "bad points" are left.
    v_points = np.zeros(void.points.shape)
    while (v_points != void.points).any() if v_points.size == void.points.size else True:
        if abs(v_points.mean().mean() - void.points.mean().mean()) < 1e-6:
            break

        v_points = np.copy(void.points)
        void.points = align_base_points(void_points=void.points)
        i_void = [(abs(void.cell_points(i)[:, 2] -
                       np.tan(alpha * np.pi / 180) * (x_max - void.cell_points(i)[:, 0])) > tol).any()
                  for i in range(void.n_cells)]
        ind_void = np.arange(void.n_cells)[i_void]
        # Extract the void.
        void = void.extract_cells(ind=ind_void)

    # Save the void surface as a .stl file.
    name = os.path.join(temp_folder, 'surf101_void.stl')
    pv.save_meshio(name, void)

    # Remesh the void surface.
    remesh_surface_gmsh(filename=name, mesh_resolution=lcv)


def find_orientation(xyz_n, xmin, xmax, ymin, ymax):
    """
    Find the indices of the points and curves relative to the ascribed bed corner points, with respect to these orders:
        Point order: [x_min, y_min] (0), [x_min, y_max] (1), [x_max, y_max] (2) and [x_max, y_min] (3).
        Curve order: P[0]-P[1] (0), P[1]-P[2] (1), P[2]-P[3] (2), P[3]-P[0] (3);
                        i.e. x=x_min (0), y=y_max (1), x=x_max (2), y=y_min (3).
        Curve orientation: "+" if P[i]-P[i+1] and "-" if P[i+i]-P[i].
    """
    # Initialize the index arrays.
    p_ind  = []     # The spatially correct order of the points.
    c_ind  = []     # The spatially correct order of the curves.
    c_sign = []     # The spatially correct orientation of the curves.
    for j in range(4):
        # Find what corner point of the top surface is directly above its relative bed counterpart.
        if xyz_n[j][0] == xmin:
            if xyz_n[j][1] == ymin:
                p_ind.append(0)
            elif xyz_n[j][1] == ymax:
                p_ind.append(1)
        elif xyz_n[j][0] == xmax:
            if xyz_n[j][1] == ymax:
                p_ind.append(2)
            elif xyz_n[j][1] == ymin:
                p_ind.append(3)

        # Find what curve is at what boundary and its current orientation(+/-).
        if xyz_n[j][0] == xyz_n[j + 1][0] if j < 3 else xyz_n[j][0] == xyz_n[0][0]:
            if xyz_n[j][0] == xmin:
                c_ind.append(0)
                c_sign.append(np.sign(xyz_n[j + 1][1] - xyz_n[j][1] if j < 3 else xyz_n[0][1] - xyz_n[j][1]))
            elif xyz_n[j][0] == xmax:
                c_ind.append(2)
                c_sign.append(
                    -1 * np.sign(xyz_n[j + 1][1] - xyz_n[j][1] if j < 3 else xyz_n[0][1] - xyz_n[j][1]))
        elif xyz_n[j][1] == xyz_n[j + 1][1] if j < 3 else xyz_n[j][1] == xyz_n[0][1]:
            if xyz_n[j][1] == ymax:
                c_ind.append(1)
                c_sign.append(np.sign(xyz_n[j + 1][0] - xyz_n[j][0] if j < 3 else xyz_n[0][0] - xyz_n[j][0]))
            elif xyz_n[j][1] == ymin:
                c_ind.append(3)
                c_sign.append(
                    -1 * np.sign(xyz_n[j + 1][0] - xyz_n[j][0] if j < 3 else xyz_n[0][0] - xyz_n[j][0]))

    return p_ind, c_ind, c_sign


def create_sides(s, c):
    """
    Create the vertical sides, as well as the bed between the walls and the void, once the void and top surface have
    been imported.
    Physical surfaces are in order: bed (101), x=x_min (102), y=y_max (103), x=x_max (104), y=y_min (105) and top (106).
    Inputs:
    s: the dimtags of the imported surfaces.
    c: the dimtags of the curves of the imported surfaces.
    """
    # Find which surface is at the top (1) and bottom (0).
    s0 = [s[i][1] for i in range(len(s)) if gmsh.model.get_closest_point(s[i][0], s[i][1], C)[0][2] < 0.75*H]
    s1 = [s[i][1] for i in range(len(s)) if gmsh.model.get_closest_point(s[i][0], s[i][1], C)[0][2] > 0.75*H]

    # Save the point numbers and their coordinates.
    p = []
    xyz = []
    for i in range(len(c)):
        pt = gmsh.model.getBoundary([c[i]], combined=False)
        p.extend([pt[0][1]])
        xyz.append(gmsh.model.getValue(0, pt[0][1], []))

    # Distinguish between points and curves on the bottom (0) and top (1).
    p0 = [p[i] for i in range(len(xyz)) if xyz[i][0] != x_min and xyz[i][0] != x_max]
    xyz0 = [xyz[i] for i in range(len(xyz)) if xyz[i][0] != x_min and xyz[i][0] != x_max]
    c0 = [c[i][1] for i in range(len(xyz)) if xyz[i][0] != x_min and xyz[i][0] != x_max]
    p1 = [p[i] for i in range(len(xyz)) if xyz[i][0] == x_min or xyz[i][0] == x_max]
    xyz1 = [xyz[i] for i in range(len(xyz)) if xyz[i][0] == x_min or xyz[i][0] == x_max]
    c1 = [c[i][1] for i in range(len(xyz)) if xyz[i][0] == x_min or xyz[i][0] == x_max]

    # Get the correct orientation of points and curves on the top surface.
    ind_p1, ind_c1, sign_c1 = find_orientation(xyz1, xmin=x_min, xmax=x_max, ymin=y_min, ymax=y_max)

    # Overwrite the arrays with the point and curve numbers in correct order (and with correct sign).
    p1 = [p1[ind_p1.index(i)] for i in range(len(p1))]
    c1 = [int(c1[ind_c1.index(i)] * sign_c1[i]) for i in range(4)]

    # Lower corner points
    gmsh.model.geo.addPoint(x_min, y_min, x_max*np.tan(alpha * np.pi / 180), lc, 11)
    gmsh.model.geo.addPoint(x_min, y_max, x_max*np.tan(alpha * np.pi / 180), lc, 12)
    gmsh.model.geo.addPoint(x_max, y_max, 0, lc, 13)
    gmsh.model.geo.addPoint(x_max, y_min, 0, lc, 14)

    # Connect the lower corner points
    gmsh.model.geo.addLine(11, 12, 1112)
    gmsh.model.geo.addLine(12, 13, 1213)
    gmsh.model.geo.addLine(13, 14, 1314)
    gmsh.model.geo.addLine(14, 11, 1411)

    # Point and curve numbers of the newly added ones (01).
    p01 = [11, 12, 13, 14]
    c01 = [1112, 1213, 1314, 1411]

    # Connect the upper and lower corners.
    for i in range(4):
        gmsh.model.geo.addLine(p01[i], p1[i], 111 * (i + 1))

    # Create curve loops for the bed.
    gmsh.model.geo.addCurveLoop(c0, 100)
    gmsh.model.geo.addCurveLoop(c01, 101)

    # Create curve loops for the vertical sides.
    gmsh.model.geo.addCurveLoop([c01[0], 222, c1[0], -111], 102)
    gmsh.model.geo.addCurveLoop([c01[1], 333, c1[1], -222], 103)
    gmsh.model.geo.addCurveLoop([c01[2], 444, c1[2], -333], 104)
    gmsh.model.geo.addCurveLoop([c01[3], 111, c1[3], -444], 105)

    # Create the bottom surface
    gmsh.model.geo.addPlaneSurface([101, 100], 101)

    # Create side surfaces.
    for i in range(4):
        gmsh.model.geo.addPlaneSurface([102 + i], 102 + i)

    # Synchronize CAD entities with the model.
    gmsh.model.geo.synchronize()

    # Add the physical entities.
    gmsh.model.geo.addPhysicalGroup(2, [101, s0[0]], 101)
    for i in range(4):
        gmsh.model.geo.addPhysicalGroup(2, [102 + i], 102 + i)
    gmsh.model.geo.addPhysicalGroup(2, s1, 106)

    # Add surface loop and volume.
    gmsh.model.geo.addSurfaceLoop([s0[0], 101, 102, 103, 104, 105, s1[0]], 1)
    gmsh.model.geo.addVolume([1], 1)
    gmsh.model.geo.addPhysicalGroup(3, [1], 1)

    # Synchronize CAD entities with the model
    gmsh.model.geo.synchronize()


def meshSizeCallback(dim, tag, x, y, z):
    """
    Compute the mesh size field. The minimum and maximum grid sizes are "lcv" and "lc", respectively.
    For points within a distance of "lc" from the void, the grid size is set equal to "lcv" but otherwise "lc".
    The function relies on the global variables: "lcv", "lc", "centre_line" and "radii".
    """
    # Read in the void mesh with pyvista.
    void = pv.read(os.path.join(temp_folder, 'surf101_void.stl'))
    # Find the smallest distance between the point and the void.
    dist = np.linalg.norm(void.points - np.array([x, y, z]), axis=1).min()

    # Return the minimum resolution if the point is within the specified distance, else the maximum resolution.
    if lc <= 15:
        if dist <= lc:
            return lcv
        else:
            return lc
    else:
        if dist <= 15:
            return lcv
        elif (dist > 15) and (dist <= 30):
            return 15
        else:
            return lc


def merge_gmsh(output_name):
    """
    Merge the void and top surface and fill in the rest of the model domain.
    """
    gmsh.initialize()
    gmsh.option.setNumber("Mesh.Algorithm", 5)
    gmsh.option.setNumber('General.Verbosity', 1)

    # Load STL files. The top surface is loaded first as it contains more boundary entities.
    gmsh.merge(os.path.join(temp_folder, 'surf106.stl'))
    gmsh.merge(os.path.join(temp_folder, 'surf101_void.stl'))
    # Classify and create geometries for the imported surfaces and their boundaries.
    gmsh.model.mesh.classifySurfaces(math.pi, curveAngle=math.pi / 4)
    gmsh.model.mesh.createGeometry()

    # retrieve the surface, its boundary curves and corner points
    s = gmsh.model.getEntities(2)
    c = gmsh.model.getBoundary(s)

    # Synchronize CAD entities with the model
    gmsh.model.geo.synchronize()

    # Fill in the rest of the domain.
    create_sides(s=s, c=c)

    # Set the mesh element size
    gmsh.model.mesh.setSizeCallback(meshSizeCallback)

    # Generate the 3D mesh.
    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.optimize('', True)

    # Save and close.
    gmsh.write(output_name)
    gmsh.finalize()


def print_message(message):
    print('  -- ' + message)


print('- Converting .vtu to .msh')
# Extract the surface meshes from the Elmer/Ice results.
print_message('Reading .vtu file and extracting meshes of all surfaces.')
vtu2stl(vtu_file=args.file)

# Extract the void from the bottom surface and remesh it.
print_message('Extracting and remeshing the void.')
extract_and_remesh_void()

# Merge the void and top surface meshes and create the 3D mesh for Elmer/Ice.
print_message('Merging surface and void and creating the mesh.')
merge_gmsh(os.path.join(cwd, geometry + '.msh'))

