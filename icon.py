import matplotlib as plt
from svgpathtools import svg2paths
from svgpath2mpl import parse_path

plane_path, attributes = svg2paths('airplane.svg')
plane_marker = parse_path(attributes[0]['d'])
plane_marker.vertices -= plane_marker.vertices.mean(axis=0)
plane_marker = plane_marker.transformed(plt.transforms.Affine2D().rotate_deg(174))
plane_marker = plane_marker.transformed(plt.transforms.Affine2D().scale(-1,1))

grass_path, attributes = svg2paths('grass.svg')
grass_marker = parse_path(attributes[0]['d'])
grass_marker.vertices -= grass_marker.vertices.mean(axis=0)
grass_marker = grass_marker.transformed(plt.transforms.Affine2D().rotate_deg(174))
grass_marker = grass_marker.transformed(plt.transforms.Affine2D().scale(-1,1))

