import gmsh
import os
import sys
import numpy as np
import argparse
import matplotlib.pyplot as plt

# Initiate the argument parser
parser = argparse.ArgumentParser()
# Add arguments
parser.add_argument('-g', '--geometry', action='store', type=str, required=True, choices=['cave', 'chimney'],
                    help='Select glaciovolcanic void geometry: "cave" or "chimney".')
parser.add_argument('-H', '--thickness', action='store', type=float, default=100,
                    help='Set glacier thickness.')
parser.add_argument('-I', '--incline', action='store', type=float, required=True,
                    help='The bed incline / the angle the gravity vector is tilted by.')
parser.add_argument('-R', '--radius', action='store', type=float,
                    help='Set the initial radial geometry of the glaciovolcanic void.')
parser.add_argument('-m', '--minRes', action='store', type=float, default=2,
                    help='Set the minimum mesh resolution around the glaciovolcanic void.')
parser.add_argument('-M', '--maxRes', action='store', type=float, default=10,
                    help='Set the maximum mesh resolution around the boundaries.')

# Read arguments from the command line
args = parser.parse_args()

# Assign the thickness based on an argument input
H = args.thickness
# Assign the slope based on an argument input.
alpha = args.incline
# Assign the initial radial geometry. A condition on the size of R is needed as it strictly has to be less than H.
if args.radius:
    R = args.radius
    if R >= H - 5:
        sys.exit('The selected radial geometry is too large! It must be less than the glacier thickness.')
else:
    R = H/10
# Assign the mesh resolution
lcv = args.minRes
lc  = args.maxRes

# Current working directory.
cwd = os.getcwd()

# Set the spatial dimensions of the domain
x_min, x_max, y_min, y_max = [0, 2*H, 0, 2*H]
C = np.array([(x_max - x_min)/2, (y_max - y_min)/2, np.tan(alpha * np.pi / 180)*(x_max - x_min)/2])


def meshSizeCallback(dim, tag, x, y, z):
    """
    Compute the mesh size field. The minimum and maximum grid sizes are "lcv" and "lc", respectively.
    For points within a distance of "max(2R,lc)" from the void, the grid size is set equal to "lcv" but otherwise "lc".
    The function relies on the global variables: "lcv" and "lc".
    """
    # Compute the distance from any given point to the centre point of the void.
    if args.geometry == "cave":
        dist = np.linalg.norm(C - np.array([x, y, z]))
    else:
        dist = np.linalg.norm(C[0:2] - np.array([x, y]))

    # Return the minimum resolution if the point is within the specified distance, else the maximum resolution.
    if lc <= 15:
        if dist <= np.array([lc, 2*R]).max():
            return lcv
        else:
            return lc
    else:
        if dist <= np.array([lc, 2*R]).max():
            return lcv
        elif (dist > np.array([lc, 2*R]).max()) and (dist <= np.array([2*lc, 4*R]).max()):
            return 15
        else:
            return lc


def deformed_cave_gmsh():
    """
    Create the model domain of a cubed glacier parcel on a sloped bed with an void prescribed at the centre of bottom.
    The void is approximated as a spherical cap of height h and radius a, cut from a sphere of the prescribed radius R.
    """
    # Start gmsh
    gmsh.initialize()
    gmsh.option.setNumber('General.Verbosity', 1)

    # Lower corner points
    gmsh.model.geo.addPoint(x_min, y_min, x_max*np.tan(alpha * np.pi / 180), lc, 11)
    gmsh.model.geo.addPoint(x_min, y_max, x_max*np.tan(alpha * np.pi / 180), lc, 12)
    gmsh.model.geo.addPoint(x_max, y_max, 0, lc, 13)
    gmsh.model.geo.addPoint(x_max, y_min, 0, lc, 14)

    # Upper corner points
    gmsh.model.geo.addPoint(x_min, y_min, H + x_max*np.tan(alpha * np.pi / 180), lc, 21)
    gmsh.model.geo.addPoint(x_min, y_max, H + x_max*np.tan(alpha * np.pi / 180), lc, 22)
    gmsh.model.geo.addPoint(x_max, y_max, H, lc, 23)
    gmsh.model.geo.addPoint(x_max, y_min, H, lc, 24)

    # Connect the corner points
    gmsh.model.geo.addLine(11, 12, 1112)
    gmsh.model.geo.addLine(12, 13, 1213)
    gmsh.model.geo.addLine(13, 14, 1314)
    gmsh.model.geo.addLine(14, 11, 1411)

    gmsh.model.geo.addLine(21, 22, 2122)
    gmsh.model.geo.addLine(22, 23, 2223)
    gmsh.model.geo.addLine(23, 24, 2324)
    gmsh.model.geo.addLine(24, 21, 2421)

    gmsh.model.geo.addLine(11, 21, 1121)
    gmsh.model.geo.addLine(12, 22, 1222)
    gmsh.model.geo.addLine(13, 23, 1323)
    gmsh.model.geo.addLine(14, 24, 1424)

    # Create the cube curve loops
    gmsh.model.geo.addCurveLoop([1112, 1213, 1314, 1411], 101)      # z = 0
    gmsh.model.geo.addCurveLoop([1112, 1222, -2122, -1121], 102)    # x = 0
    gmsh.model.geo.addCurveLoop([1213, 1323, -2223, -1222], 103)    # y = y_max
    gmsh.model.geo.addCurveLoop([1314, 1424, -2324, -1323], 104)    # x = x_max
    gmsh.model.geo.addCurveLoop([1411, 1121, -2421, -1424], 105)    # y = 0
    gmsh.model.geo.addCurveLoop([2122, 2223, 2324, 2421], 106)      # z = H

    # Create the cave, assumed here to be a spherical cap of radius a and height h cut from a sphere of radius R.
    cos = np.cos(alpha * np.pi / 180)
    sin = np.sin(alpha * np.pi / 180)
    h = R / 2                       # The height of the spherical cap.
    a = np.sqrt(R**2 - (R-h)**2)    # The radius of the spherical cap.
    gmsh.model.geo.addPoint(C[0], C[1], C[2], lcv, 30)                          # Centre
    gmsh.model.geo.addPoint(C[0], C[1] - a, C[2], lcv, 31)                      # East
    gmsh.model.geo.addPoint(C[0] - cos*a, C[1], C[2] + sin*a, lcv, 32)          # South
    gmsh.model.geo.addPoint(C[0], C[1] + a, C[2], lcv, 33)                      # West
    gmsh.model.geo.addPoint(C[0] + cos*a, C[1], C[2] - sin*a, lcv, 34)          # North
    gmsh.model.geo.addPoint(C[0] + sin*h, C[1], C[2] + cos*h, lcv, 35)          # Top
    gmsh.model.geo.addPoint(C[0] - sin*(R-h), C[1], C[2] - cos*(R-h), lcv, 36)  # Bottom of sphere

    # Connect the cave points (into quarter-length circle arcs)
    gmsh.model.geo.addCircleArc(31, 30, 32, 3132)  # E-S
    gmsh.model.geo.addCircleArc(32, 30, 33, 3233)  # S-W
    gmsh.model.geo.addCircleArc(33, 30, 34, 3334)  # W-N
    gmsh.model.geo.addCircleArc(34, 30, 31, 3431)  # N-E

    gmsh.model.geo.addCircleArc(31, 36, 35, 3135)  # E-T
    gmsh.model.geo.addCircleArc(32, 36, 35, 3235)  # S-T
    gmsh.model.geo.addCircleArc(33, 36, 35, 3335)  # W-T
    gmsh.model.geo.addCircleArc(34, 36, 35, 3435)  # N-T

    # Connect the cave arcs (hemisphere quarters)
    gmsh.model.geo.addCurveLoop([3132, 3233, 3334, 3431], 300)

    gmsh.model.geo.addCurveLoop([3132, 3235, -3135], 301)
    gmsh.model.geo.addCurveLoop([3233, 3335, -3235], 302)
    gmsh.model.geo.addCurveLoop([3334, 3435, -3335], 303)
    gmsh.model.geo.addCurveLoop([3431, 3135, -3435], 304)

    # Create the surfaces
    gmsh.model.geo.addPlaneSurface([101, 300], 101)
    gmsh.model.geo.addPlaneSurface([102], 102)
    gmsh.model.geo.addPlaneSurface([103], 103)
    gmsh.model.geo.addPlaneSurface([104], 104)
    gmsh.model.geo.addPlaneSurface([105], 105)
    gmsh.model.geo.addPlaneSurface([106], 106)

    gmsh.model.geo.addSurfaceFilling([301], 301)
    gmsh.model.geo.addSurfaceFilling([302], 302)
    gmsh.model.geo.addSurfaceFilling([303], 303)
    gmsh.model.geo.addSurfaceFilling([304], 304)

    # Create the volume
    gmsh.model.geo.addSurfaceLoop([301, 302, 303, 304, 101, 102, 103, 104, 105, 106], 1)
    gmsh.model.geo.addVolume([1], 1)

    # Synchronize CAD entities with the model
    gmsh.model.geo.synchronize()

    # Create the physical entities
    gmsh.model.addPhysicalGroup(2, [101, 301, 302, 303, 304], 101)
    gmsh.model.addPhysicalGroup(2, [102], 102)
    gmsh.model.addPhysicalGroup(2, [103], 103)
    gmsh.model.addPhysicalGroup(2, [104], 104)
    gmsh.model.addPhysicalGroup(2, [105], 105)
    gmsh.model.addPhysicalGroup(2, [106], 106)

    gmsh.model.addPhysicalGroup(3, [1], 1)

    # Set the mesh element size
    gmsh.model.mesh.setSizeCallback(meshSizeCallback)

    # Generate the 3D mesh
    gmsh.model.mesh.generate(3)

    # Optimize the mesh
    gmsh.model.mesh.optimize('', True)

    # Export the mesh
    gmsh.write(os.path.join(cwd, 'cave.msh'))
    gmsh.finalize()


def deformed_chimney_gmsh():
    """
    Create the model domain of a cubed glacier parcel on a sloped bed with a void prescribed from the centre of bottom.
    The void is approximated as cylinder, of radius R, topped by a hemisphere.
        The chimney stretches up to 5 m below the surface.
    """
    # Start gmsh
    gmsh.initialize()
    gmsh.option.setNumber('General.Verbosity', 1)

    # Lower corner points
    gmsh.model.geo.addPoint(x_min, y_min, x_max*np.tan(alpha * np.pi / 180), lc, 11)
    gmsh.model.geo.addPoint(x_min, y_max, x_max*np.tan(alpha * np.pi / 180), lc, 12)
    gmsh.model.geo.addPoint(x_max, y_max, 0, lc, 13)
    gmsh.model.geo.addPoint(x_max, y_min, 0, lc, 14)

    # Upper corner points
    gmsh.model.geo.addPoint(x_min, y_min, H + x_max*np.tan(alpha * np.pi / 180), lc, 21)
    gmsh.model.geo.addPoint(x_min, y_max, H + x_max*np.tan(alpha * np.pi / 180), lc, 22)
    gmsh.model.geo.addPoint(x_max, y_max, H, lc, 23)
    gmsh.model.geo.addPoint(x_max, y_min, H, lc, 24)

    # Connect the corner points
    gmsh.model.geo.addLine(11, 12, 1112)
    gmsh.model.geo.addLine(12, 13, 1213)
    gmsh.model.geo.addLine(13, 14, 1314)
    gmsh.model.geo.addLine(14, 11, 1411)

    gmsh.model.geo.addLine(21, 22, 2122)
    gmsh.model.geo.addLine(22, 23, 2223)
    gmsh.model.geo.addLine(23, 24, 2324)
    gmsh.model.geo.addLine(24, 21, 2421)

    gmsh.model.geo.addLine(11, 21, 1121)
    gmsh.model.geo.addLine(12, 22, 1222)
    gmsh.model.geo.addLine(13, 23, 1323)
    gmsh.model.geo.addLine(14, 24, 1424)

    # Create the cube curve loops
    gmsh.model.geo.addCurveLoop([1112, 1213, 1314, 1411], 101)      # z = 0
    gmsh.model.geo.addCurveLoop([1112, 1222, -2122, -1121], 102)    # x = 0
    gmsh.model.geo.addCurveLoop([1213, 1323, -2223, -1222], 103)    # y = y_max
    gmsh.model.geo.addCurveLoop([1314, 1424, -2324, -1323], 104)    # x = x_max
    gmsh.model.geo.addCurveLoop([1411, 1121, -2421, -1424], 105)    # y = 0
    gmsh.model.geo.addCurveLoop([2122, 2223, 2324, 2421], 106)      # z = H

    # Create the chimney, assumed here to be cylinder of radius R, topped by a hemisphere.
    z_C = np.tan(alpha * np.pi / 180) * R  # The difference in height between centre and N/S.
    dh = lcv                               # The cap thickness between void and top surface.
    h = H - R - dh                         # The height of chimney stem.

    # The chimney points at the bed.
    gmsh.model.geo.addPoint(C[0], C[1], C[2], lcv, 30)            # Centre
    gmsh.model.geo.addPoint(C[0], C[1] - R, C[2], lcv, 31)        # East
    gmsh.model.geo.addPoint(C[0] - R, C[1], C[2] + z_C, lcv, 32)  # South
    gmsh.model.geo.addPoint(C[0], C[1] + R, C[2], lcv, 33)        # West
    gmsh.model.geo.addPoint(C[0] + R, C[1], C[2] - z_C, lcv, 34)  # North

    # The chimney points at the top.
    gmsh.model.geo.addPoint(C[0], C[1], C[2] + h, lcv, 40)      # Centre
    gmsh.model.geo.addPoint(C[0], C[1] - R, C[2] + h, lcv, 41)  # East
    gmsh.model.geo.addPoint(C[0] - R, C[1], C[2] + h, lcv, 42)  # South
    gmsh.model.geo.addPoint(C[0], C[1] + R, C[2] + h, lcv, 43)  # West
    gmsh.model.geo.addPoint(C[0] + R, C[1], C[2] + h, lcv, 44)  # North
    gmsh.model.geo.addPoint(C[0], C[1], C[2] + h + R, lcv, 45)  # Top

    # Connect the chimney base points (into quarter-length ellipse arcs)
    gmsh.model.geo.addEllipseArc(31, 30, 32, 32, 3132)  # E-S
    gmsh.model.geo.addEllipseArc(32, 30, 32, 33, 3233)  # S-W
    gmsh.model.geo.addEllipseArc(33, 30, 34, 34, 3334)  # W-N
    gmsh.model.geo.addEllipseArc(34, 30, 34, 31, 3431)  # N-E

    # Connect the chimney top points (into quarter-length circle arcs)
    gmsh.model.geo.addCircleArc(41, 40, 42, 4142)  # E-S
    gmsh.model.geo.addCircleArc(42, 40, 43, 4243)  # S-W
    gmsh.model.geo.addCircleArc(43, 40, 44, 4344)  # W-N
    gmsh.model.geo.addCircleArc(44, 40, 41, 4441)  # N-E

    gmsh.model.geo.addCircleArc(41, 40, 45, 4145)  # E-T
    gmsh.model.geo.addCircleArc(42, 40, 45, 4245)  # S-T
    gmsh.model.geo.addCircleArc(43, 40, 45, 4345)  # W-T
    gmsh.model.geo.addCircleArc(44, 40, 45, 4445)  # N-T

    # Connect the base and top of chimney.
    gmsh.model.geo.addLine(31, 41, 3141)
    gmsh.model.geo.addLine(32, 42, 3242)
    gmsh.model.geo.addLine(33, 43, 3343)
    gmsh.model.geo.addLine(34, 44, 3444)

    # Connect the base arcs (hemisphere quarters)
    gmsh.model.geo.addCurveLoop([3132, 3233, 3334, 3431], 300)
    # Connect the stem sectors.
    gmsh.model.geo.addCurveLoop([3132, 3242, -4142, -3141], 301)
    gmsh.model.geo.addCurveLoop([3233, 3343, -4243, -3242], 302)
    gmsh.model.geo.addCurveLoop([3334, 3444, -4344, -3343], 303)
    gmsh.model.geo.addCurveLoop([3431, 3141, -4441, -3444], 304)
    # Connect the top sectors.
    gmsh.model.geo.addCurveLoop([4142, 4245, -4145], 401)
    gmsh.model.geo.addCurveLoop([4243, 4345, -4245], 402)
    gmsh.model.geo.addCurveLoop([4344, 4445, -4345], 403)
    gmsh.model.geo.addCurveLoop([4441, 4145, -4445], 404)

    # Create the surfaces
    gmsh.model.geo.addPlaneSurface([101, 300], 101)
    gmsh.model.geo.addPlaneSurface([102], 102)
    gmsh.model.geo.addPlaneSurface([103], 103)
    gmsh.model.geo.addPlaneSurface([104], 104)
    gmsh.model.geo.addPlaneSurface([105], 105)
    gmsh.model.geo.addPlaneSurface([106], 106)

    gmsh.model.geo.addSurfaceFilling([301], 301)
    gmsh.model.geo.addSurfaceFilling([302], 302)
    gmsh.model.geo.addSurfaceFilling([303], 303)
    gmsh.model.geo.addSurfaceFilling([304], 304)

    gmsh.model.geo.addSurfaceFilling([401], 401)
    gmsh.model.geo.addSurfaceFilling([402], 402)
    gmsh.model.geo.addSurfaceFilling([403], 403)
    gmsh.model.geo.addSurfaceFilling([404], 404)

    # Create the volume
    gmsh.model.geo.addSurfaceLoop([401, 402, 403, 404, 304, 303, 302, 301, 101, 102, 103, 104, 105, 106], 1)
    gmsh.model.geo.addVolume([1], 1)

    # Synchronize CAD entities with the model
    gmsh.model.geo.synchronize()

    # Create the physical entities
    gmsh.model.addPhysicalGroup(2, [101, 301, 302, 303, 304, 404, 403, 402, 401], 101)
    gmsh.model.addPhysicalGroup(2, [102], 102)
    gmsh.model.addPhysicalGroup(2, [103], 103)
    gmsh.model.addPhysicalGroup(2, [104], 104)
    gmsh.model.addPhysicalGroup(2, [105], 105)
    gmsh.model.addPhysicalGroup(2, [106], 106)

    gmsh.model.addPhysicalGroup(3, [1], 1)

    # Set the mesh element size
    gmsh.model.mesh.setSizeCallback(meshSizeCallback)

    # Generate the 3D mesh
    gmsh.model.mesh.generate(3)

    # Optimize the mesh
    gmsh.model.mesh.optimize('', True)

    # Export the mesh
    gmsh.write(os.path.join(cwd, 'chimney.msh'))
    gmsh.finalize()


def datafile(x0, xe, y0, ye):
    # Set the resolution to 0.5 m and find the number of points needed.
    nx = int((xe - x0 + 2)/0.5 + 1)
    ny = int((ye - y0 + 2)/0.5 + 1)

    # Create the arrays for the bed grid.
    x = np.linspace(x0-1, xe+1, nx)
    y = np.linspace(y0-1, ye+1, ny)

    # Grid the data.
    X, Y = np.meshgrid(x, y)
    Z = np.tan(alpha * np.pi / 180) * (xe - X)

    # Only save 3 decimals
    X = np.round(X, 3)
    Y = np.round(Y, 3)
    Z = np.round(Z, 3)

    # Create/open the file to write
    data = open('Z_bed.dat', 'w')

    # Write into the file.
    for i in range(0, len(X[:, 0])):
        for j in range(0, len(X[0, :])):
            data.write('{0} {1} {2}\n'.format(X[i, j], Y[i, j], Z[i, j]))


# Create the gmsh files
if args.geometry == 'cave':
    deformed_cave_gmsh()
elif args.geometry == 'chimney':
    deformed_chimney_gmsh()

# Create a .dat file to save a gridded representation of the bed.
datafile(x_min, x_max, y_min, y_max)
