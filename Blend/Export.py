import bpy
import json
import mathutils
import os

def get_camera_data(camera):
    """Get camera data"""
    location = camera.location
    camera_matrix = camera.matrix_world
    forward = (camera_matrix.to_3x3() @ mathutils.Vector((0, 0, -1))).normalized()
    gaze_vector = {
        "x": forward.x,
        "y": forward.y,
        "z": forward.z
    }
    up_vector_world = (camera_matrix.to_3x3() @ mathutils.Vector((0, 1, 0))).normalized()
    up_vector = {
        "x": up_vector_world.x,
        "y": up_vector_world.y,
        "z": up_vector_world.z
    }

    aperture = camera.data.dof.aperture_fstop
    focal_length = camera.data.lens
    focal_distance = getattr(camera.data.dof, 'focus_distance', None)
    
    sensor = {
        "width": camera.data.sensor_width,
        "height": camera.data.sensor_height
    }
    
    scene = bpy.context.scene
    
    film_resolution = {
        "width": scene.render.resolution_x,
        "height": scene.render.resolution_y
    }
    
    return {
        "location": {"x": location.x, "y": location.y, "z": location.z},
        "gaze_vector": gaze_vector,
        "up_vector": up_vector,
        "aperture": aperture,
        "focal_distance": focal_distance,
        "focal_length": focal_length,
        "sensor": sensor,
        "film_resolution": film_resolution
    }



def get_point_light_data(point_light):
    """Get point light data"""
    location = point_light.location
    radiant_intensity = point_light.data.energy
    return {
        "location": {"x": location.x, "y": location.y, "z": location.z},
        "radiant_intensity": radiant_intensity
    }

def get_sphere_data(sphere):
    """Get sphere data"""
    scene = bpy.context.scene
    start_frame = scene.frame_start
    end_frame = scene.frame_end
    
    start_location = sample_object_transform(sphere, start_frame)
    end_location = sample_object_transform(sphere, end_frame)
    
    radius = sphere.scale.x  # spheres are often scaled equally in all directions
    material = get_material_data(sphere)
    
    return {
        "start_location": {"x": start_location["x"], "y": start_location["y"], "z": start_location["z"]},
        "end_location": {"x": end_location["x"], "y": end_location["y"], "z": end_location["z"]},
        "radius": radius,
        "material": material
    }

def get_cube_data(cube):
    """Get cube data"""
    scene = bpy.context.scene
    start_frame = scene.frame_start
    end_frame = scene.frame_end
    
    start_location = sample_object_transform(cube, start_frame)
    end_location = sample_object_transform(cube, end_frame)
    
    rotation = cube.rotation_euler
    scale = cube.scale.x  # cubes are often scaled equally in all directions
    material = get_material_data(cube)
    
    return {
        "start_location": {"x": start_location["x"], "y": start_location["y"], "z": start_location["z"]},
        "end_location": {"x": end_location["x"], "y": end_location["y"], "z": end_location["z"]},
        "rotation": {"x": rotation.x, "y": rotation.y, "z": rotation.z},
        "scale": scale,
        "material": material
    }

def get_plane_data(plane):
    """Get plane data"""
    # Get the corners of the plane by applying the matrix_world to the vertices
    corners = []
    material = get_material_data(plane)
    for vertex in plane.data.vertices:
        corner = plane.matrix_world @ vertex.co
        corners.append({"x": corner.x, "y": corner.y, "z": corner.z})
    return {
        "corners": corners,
        "material": material
    }


def get_material_data(obj):
    """Extract material info from a mesh object"""
    if not obj.data.materials or obj.data.materials[0] is None:
        # Default material if none assigned
        return {
            "diffuse": {"r": 1.0, "g": 1.0, "b": 1.0},
            "specular": {"r": 0.3, "g": 0.3, "b": 0.3},
            "shininess": 32.0,
            "transparency": 0.0,
            "ior": 1.0,
            "texture": None
        }

    mat = obj.data.materials[0]  # Take the first material for simplicity

    # Default values
    diffuse = {"r": 1.0, "g": 1.0, "b": 1.0}
    specular = {"r": 0.3, "g": 0.3, "b": 0.3}
    shininess = 32.0
    transparency = 0.0
    ior = 1.0
    texture_path = None
    
    principled_found = False

    if mat.node_tree:
        for node in mat.node_tree.nodes:
            if node.type == 'BSDF_PRINCIPLED':
                principled_found = True

                color_input = node.inputs['Base Color']
                base_color = color_input.default_value
                diffuse = {"r": base_color[0], "g": base_color[1], "b": base_color[2]}

                # Traverse texture if linked
                if color_input.is_linked:
                    from_node = color_input.links[0].from_node
                    tex = find_connected_texture(from_node)
                    if tex:
                        texture_path = tex

                # Specular
                if 'Specular' in node.inputs:
                    s = node.inputs['Specular'].default_value
                    s_clamped = min(s, 0.3)
                    specular = {"r": s_clamped, "g": s_clamped, "b": s_clamped}

                # Roughness → shininess
                if 'Roughness' in node.inputs:
                    roughness = node.inputs['Roughness'].default_value
                    shininess = max(0.0, (1.0 - roughness) * 128.0)

                # Transparency
                if 'Transmission' in node.inputs:
                    transmission = node.inputs['Transmission'].default_value
                    transparency = max(0.0, min(1.0, transmission))

                # Refractive index
                if 'IOR' in node.inputs:
                    ior = node.inputs['IOR'].default_value


                break  # Only process first Principled BSDF node
            
        if not principled_found:
            for node in mat.node_tree.nodes:
                if node.type == 'BSDF_DIFFUSE':
                    if 'Color' in node.inputs:
                        color_node = node.inputs['Color']
                        if hasattr(color_node, "default_value"):
                            base_color = color_node.default_value
                            diffuse = {"r": base_color[0], "g": base_color[1], "b": base_color[2]}
                            specular = {"r": 0.5, "g": 0.5, "b": 0.5}
                            shininess = 32.0
                            break
            

    return {
        "diffuse": diffuse,
        "specular": specular,
        "shininess": shininess,
        "transparency": transparency,
        "ior": ior,
        "texture": texture_path if texture_path else None
    }


def get_material_texture(mat):
    """Return the filename of the first image texture in the material, if any."""
    if not mat.node_tree:
        return None

    for node in mat.node_tree.nodes:
        tex = find_connected_texture(node)
        if tex:
            return tex

    return None  # No texture found

def find_connected_texture(node, visited=None):
    """Recursively traverse all node inputs to find the first image texture."""
    if visited is None:
        visited = set()
    if node in visited:
        return None
    visited.add(node)

    if node.type == 'TEX_IMAGE' and node.image:
        return os.path.basename(node.image.filepath_raw) or node.image.name

    if hasattr(node, 'inputs'):
        for input in node.inputs:
            if input.is_linked:
                for link in input.links:
                    from_node = link.from_node
                    tex = find_connected_texture(from_node, visited)
                    if tex:
                        return tex
    return None




def sample_object_transform(obj, frame):
    scene = bpy.context.scene
    scene.frame_set(frame)

    loc = obj.matrix_world.to_translation()

    return {
        "x": loc.x, "y": loc.y, "z": loc.z,
    }


def export_to_json():
    """Export the scene data to JSON"""
    
    
    start = bpy.context.scene.frame_start
    end = bpy.context.scene.frame_end
    
    bpy.context.scene.frame_set(start)
    
    cameras = []
    point_lights = []
    spheres = []
    cubes = []
    planes = []

    for object in bpy.context.scene.objects:
        
        print(object.data.name)
        
        if object.type == 'CAMERA':
            cameras.append(get_camera_data(object))
        elif object.type == 'LIGHT' and object.data.type == 'POINT':
            point_lights.append(get_point_light_data(object))
        elif (object.type == 'MESH' and object.data.name.startswith('Sphere')) or object.type == "META":
            spheres.append(get_sphere_data(object))
        elif object.type == 'MESH' and object.data.name.startswith('Cube'):
            cubes.append(get_cube_data(object))
        elif object.type == 'MESH' and object.data.name.startswith('Plane'):
            planes.append(get_plane_data(object))

    data = {
        "properties": {
            "cameras": cameras,
            "point_lights": point_lights,
            "spheres": spheres,
            "cubes": cubes,
            "planes": planes
        }
    }
    
    # Get the path to the .blend file
    blend_file_path = bpy.data.filepath
    # Get the directory of the .blend file
    blend_file_dir = os.path.dirname(blend_file_path)

    # Export to JSON
    with open(os.path.join(blend_file_dir, '..\\ASCII\scene.json'), 'w') as f:
        json.dump(data, f, indent=4)
        
    bpy.context.scene.frame_set(start)

export_to_json()
