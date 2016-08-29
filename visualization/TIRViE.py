#
# ---Prerequisites note---
#
# To use this you need to have the following packages installed:
#
# - VPython, usually distributed as linux package 'python-visual'
# - ImageMagick, usually contained by default in modern Linux distributions
#
# On most Linux systems the following should do the job of installing them:
#
#   $ sudo apt-get install python-visual imagemagick
#
#
import sys
import os
import random
import time
import numpy as np

import visual as v
import colorsys # only needed for drawing tomeks_double_helix() objects

from _gfrd import Sphere, Cylinder, Box, Plane, Disk
from multi import Multi
from single import *
from pair import *

__all__ = [
    'TIRViE',
    'TheImperialRoyalVisualizerForEGFRD',    
    ]

    

# Define a small number that will serve to define "zero extensions"
ALMOST_ZERO = 1e-20
    


class TheImperialRoyalVisualizerForEGFRD:
  """

    ---The Imperial-Royal Visualizer for eGFRD---
    
    started upon urgent request by his Majesty the Emperor
    by Geh. Computationsrat Dr. T.R. Sokolowski (tomek@hof.kuk)
    ;-)
    
    Version 0.95
    
    Wien, 2016
        
    TODO Add more info about usage etc.
    
  """
  def __init__(self, eGFRDsim, particle_color_dict={}, structure_color_dict={}, capture=False, capture_dir='./tirvie_frames'):

    self.sim = eGFRDsim
    self.world_size = self.get_world_size()
    
    self.objects = []

    self.scene_width  = 1000
    self.scene_height = 1000
    self.window_title = "TIRViE -- The Imperial-Royal Visualizer for eGFRD"
    self.window_id    = '' # set later by self.get_window_id(), once a scene has been created
    
    self.SCENE_BACKGROUND_COLOR = v.color.white
    
    self.info_label_strings_order = ('b', 's', 'd', 'r', 'c', 'o', 'w', 'i')
    self.INFO_TIME_ONLY           = False
    
    # Some flags tuning visual behavior of the objects
    self.SHOW_SIM_BOX     = True
    self.WIREFRAME_BOX    = False
    self.WIREFRAME_PLANES = False
    self.SHOW_DOMAINS     = True
    self.SHOW_STRUCTURES  = True
    self.RODS_AS_HELICES  = False
    self.SIM_BOX_OPACITY  = 0.05 # sth. around 0.01 -- 0.10 works well
    self.SHELL_OPACITY    = 0.10 
    self.PLANAR_SURFACES_OPACITY   = 0.20
    self.CYLINDER_SURFACES_OPACITY = 0.20
    self.DISK_SURFACES_OPACITY     = 0.60
    self.PARTICLE_RADIUS_FACTOR    = 1  # will show particle sphere radii this times larger than actual
    self.ROD_RADIUS_FACTOR         = 2  # will show radii of cyl. surfaces this times larger than actual
    self.N_HELIX_COILS             = 30 # sth. btw. 15-30 seems to be a good choice; this is eyecandy anyhow ;-)
    self.WIREFRAME_RADIUS          = self.world_size/500 # seems fair enough
    self.HELIX_THICKNESS_FACTOR    = 0.4
    
    # No. of random colors in addition to primary
    self.N_RANDOM_COLORS = 1000
    
    # A flag regulating output capturing
    self.CAPTURE = capture
    # The directory name where the captured frames will be stored
    self.capture_dir = capture_dir
    self.output_image_suffix = '.png' # should be in format .xxx (.png, .jpg, etc.)
    
    # A flag to pause/resume the execution of the simulation
    self.PAUSED = False
    
    # A flag to activate/inactivate retarding of progress in case
    # the visualization cannot keep speed with the simulation
    self.RETARD = False
    self.RETARDING_SECONDS = 0.1
    
    # Color lists and dictionaries for objects of different species    
    self.create_default_color_list()
    
    if particle_color_dict == {}:
      self.create_particle_color_dict()
    else:
      self.particle_color_dict = particle_color_dict
    
    if structure_color_dict == {}:
      self.create_structure_color_dict()
    else:
      self.structure_color_dict = structure_color_dict
    
    self.create_shell_color_dict()

    # Prepare the scene
    self.create_new_scene()
    
    # Create output capturing directory
    self.create_output_directory()      
        
        
  #################################################################        
  ### TESTING #####################################################
  #################################################################
  def visualize(self):    
    
    # Ask for pressed keys and handle them
    if self.scene.kb.keys:
      key = self.scene.kb.getkey()
      self.handle_key(key)

    # Ask for mouse click and handle it
    if self.scene.mouse.events:
      m = self.scene.mouse.getevent()      
      if m.click:
        # Move the info label to the position where the user clicked
        self.info_label.pos = m.pos
        
    if self.PAUSED:
      self.pause()

    if self.RETARD:
      self.retard()
      
    self.clear_scene()    
    
    Bp, Bo, Bs, Bc = self.get_cuboidal_region_data()
    self.draw_cuboidal_region(Bp, Bo, Bs)
    
    Pp, Po, Ps, Pc = self.get_planar_surface_data()
    self.draw_planar_surfaces(Pp, Po, Ps, Pc)
    
    Cp, Co, Cr, Cl, Cc = self.get_cylindrical_surface_data()
    self.draw_cylindrical_surfaces(Cp, Co, Cr, Cl, Cc)
    
    Pp, Pr, Pc = self.get_particle_data()
    self.draw_particles(Pp, Pr, Pc)
    
    sph_shell_data, cyl_shell_data = self.get_protective_domain_data()    
    self.draw_protective_domains(sph_shell_data, cyl_shell_data)        

    self.update_info_label()
    
    if self.CAPTURE:
    
      self.capture()
      
    
    #sys.exit()
  #################################################################
  #################################################################
  #################################################################
  
  
  ####################################      
  ### COLOR LISTS AND DICTIONARIES ###
  ####################################
  def create_default_color_list(self):
    
    self.default_color_list = []
    
    # First colors are the primary ones
    for c in [v.color.red, v.color.green, v.color.blue, v.color.cyan, v.color.magenta, v.color.yellow, v.color.orange]:
      
      self.default_color_list.append(c)
      
    # From there we create random RGB colors
    for n in range(self.N_RANDOM_COLORS):
        R = random.random()
        G = random.random()
        B = random.random()
        
        self.default_color_list.append( (R,G,B) )      
      
    
  def create_particle_color_dict(self):   
    
    self.particle_color_dict = {}
    
    n = 0
    for species in self.sim.world.species:
      
          c = self.default_color_list[n]
          key = species.id.serial
          self.particle_color_dict[key] = c
          n = n + 1

  
  def create_structure_color_dict(self):   
    
    self.structure_color_dict = {}
    
    n = 0
    for struct_type in self.sim.world.structure_types:
      
          c = self.default_color_list[n]
          key = struct_type.id.serial
          self.structure_color_dict[key] = c
          n = n + 1

          
  def create_shell_color_dict(self):
    
    self.shell_color_dict = {}
    
    # The "top level" shell, i.e. the one to be updated next
    self.shell_color_dict[0] = v.color.magenta
    
    # 'Single' shells
    self.shell_color_dict[1] = v.color.blue
    
    # 'Pair' shells
    self.shell_color_dict[2] = v.color.green
    
    # 'Multi" shells
    self.shell_color_dict[3] = v.color.orange
    
  
  ######################################      
  ### SCENE / OUTPUT WINDOW HANDLING ###
  ######################################
  def create_new_scene(self):
    
    try:
      #self.scene.delete()
      window.delete_all()
      self.objects = []
    except:
      print "No scene existing yet, creating now."
    
    self.scene = v.display(width = self.scene_width, height = self.scene_height, title = self.window_title)
    self.scene.background = self.SCENE_BACKGROUND_COLOR
    self.scene.select()
    
    self.info_label   = v.label(pos=(-1.1*self.world_size, 1.1*self.world_size, 0), text = 't = 0', height=10)
    self.paused_label = v.label(text = 'PAUSED', color = v.color.red, linecolor = v.color.red, visible = False)
    
    #self.logo = v.text(pos=(-self.world_size, 1.1*self.world_size, 0), text='TIRViE', height=0.1*self.world_size, width=0.2*self.world_size)
    # Sadly not working with older versions of vpython...
    
    self.get_window_id()
    
    print "Scene created."    
    
  
  def clear_scene(self, types_list = [v.box, v.sphere, v.cylinder, v.ring, v.helix, v.frame, v.curve]):
      
    for t in types_list:
      for obj in self.objects:
          if isinstance(obj, t):
            
            # First make the object invisible            
            obj.visible = False
            
            # Now delete the object
            del obj
        
    self.objects = []
    

  #########################      
  ### CAPTURING METHODS ###
  #########################
  def get_window_id(self):
    
    shell_command = 'xprop -name "' + self.window_title + '" | grep "COLORMAP" | cut -d \# -f 2 | tr -d " "'
    self.window_id = str( os.popen(shell_command).read() )    
    
    
  def create_output_directory(self):
    
    try:
      os.mkdir(self.capture_dir)
    except OSError:
      print "WARNING: Cannot create (already existant?) TIRViE output directory."
      
    
  def capture(self):
    
    # Check how many output frames there are already
    frame_no = int( os.popen('ls ' + self.capture_dir + ' | wc -l').read() )
    # Increase the frame_no by one and give it the right string format
    formated_frame_no = "%06u" % (frame_no + 1)    
    # Give the new frame the corresponding name
    savename = self.capture_dir + '/frame_' + formated_frame_no + self.output_image_suffix
    # Grab the image (using ImageMagicx) and save it
    print 'Writing ', savename
    os.popen('import -window "' + self.window_title + '" ' + savename)        
    
    
  ####################################      
  ### KEY HANDLING AND INFO OUTPUT ###
  #################################### 
  def handle_key(self, k):  
    
    # Pausing
    if k == 'p' or k == 'P':
      
      if self.PAUSED == False:
        
        self.PAUSED = True
        self.paused_label.visible = True
        print '***** PAUSED *****'
        
      else:
      
        self.PAUSED = False
        self.paused_label.visible = False
        print '***** CONTINUING *****'

    # Screen output retarding
    elif k == 'r' or k == 'R':
      
      self.RETARD = not(self.RETARD)
      
    # Visibility of info label
    if k == 'i' or k == 'I':
      
      self.info_label.visible = not(self.info_label.visible)
      
    # Show only time in info label
    if k == 't' or k == 'T':
      
      self.INFO_TIME_ONLY = not(self.INFO_TIME_ONLY)
        
    # Visibility of simulation box
    elif k == 'b' or k == 'B':
      
      if self.SHOW_SIM_BOX == False:
        
        self.SHOW_SIM_BOX = True
        print '***** SHOWING SIMULATION BOX ("world") *****'
        
      else:
        
        self.SHOW_SIM_BOX = False
        print '***** HIDING SIMULATION BOX ("world") *****'

    # Visibility of protective domains
    elif k == 'd' or k == 'D':
      
      if self.SHOW_DOMAINS == False:
        
        self.SHOW_DOMAINS = True
        print '***** SHOWING PROTECTIVE DOMAIN SHELLS *****'
        
      else:
        
        self.SHOW_DOMAINS = False
        print '***** HIDING PROTECTIVE DOMAIN SHELLS *****'
        
    # Visibility of static structures
    elif k == 's' or k == 'S':
      
      if self.SHOW_STRUCTURES == False:
        
        self.SHOW_STRUCTURES = True
        print '***** SHOWING STATIC STRUCTURES *****'
        
      else:
      
        self.SHOW_STRUCTURES = False
        print '***** HIDING STATIC STRUCTURES *****'

    # Style of drawing the cylindrical structures / rods
    elif k == 'h' or k == 'H':
      
      if self.RODS_AS_HELICES == False:
        
        self.RODS_AS_HELICES = True
        print '***** DRAWING RODS AS HELICES *****'
        
      else:
      
        self.RODS_AS_HELICES = False
        print '***** DRAWING RODS AS CYLINDERS *****'
        
    # Style of drawing the simulation box
    elif k == 'w' or k == 'W':
      
      self.WIREFRAME_BOX = not(self.WIREFRAME_BOX)
      
    # Style of drawing the planar surfaces
    elif k == 'q' or k == 'Q':
      
      self.WIREFRAME_PLANES = not(self.WIREFRAME_PLANES)
        
    # Toggling output capturing
    elif k == 'c' or k == 'C':
      
      self.CAPTURE = not(self.CAPTURE)


  def pause(self):
  
    while self.PAUSED:
      
      time.sleep(0.25)
        # prevents permanent key callback query and thus allows for better scrolling
      
      if self.scene.kb.keys:
        key = self.scene.kb.getkey()
        self.handle_key(key)


  def retard(self):
    
      time.sleep(self.RETARDING_SECONDS)
      
      
  def set_info_label_strings(self):
    
    self.info_label_strings = {}        
      
    if self.SHOW_SIM_BOX:
      self.info_label_strings['b'] = '\n[b] showing sim. box'
    else:
      self.info_label_strings['b'] = '\n[b] hiding sim. box'      
      
    if self.SHOW_DOMAINS:
      self.info_label_strings['d'] = '\n[d] showing domain shells'
    else:
      self.info_label_strings['d'] = '\n[d] hiding domain shells'
      
    if self.SHOW_STRUCTURES:
      self.info_label_strings['s'] = '\n[s] showing structures'
    else:
      self.info_label_strings['s'] = '\n[s] hiding structures'
      
    if self.RETARD:
      self.info_label_strings['r'] = '\n[r] output retarding ON'
    else:
      self.info_label_strings['r'] = '\n[r] output retarding OFF'
      
    if self.CAPTURE:
      self.info_label_strings['c'] = '\n[c] output capturing ON'
    else:
      self.info_label_strings['c'] = '\n[c] output capturing OFF'                        
    
    #if self.RODS_AS_HELICES:
      #self.info_label_strings['h'] = '\nrods shown as helices'
    #else:
      #self.info_label_strings['h'] = '\nrods shown as cylinders'
    
    self.info_label_strings['o'] = '\n\n' + str(len(self.scene.objects)) + ' objects on scene'
    
    self.info_label_strings['w'] = '\n\npress [w]/[q] for wireframe box / plane'
    self.info_label_strings['i'] = '\npress [i]/[t] to (un)hide / show time only\nleft-click to move info label position'
    
      
  def update_info_label(self):
    
    self.set_info_label_strings()
    
    try:
      time_str = "%0.3f" % self.sim.t
      newtext = 't = ' + time_str
      
      if not(self.INFO_TIME_ONLY):
        newtext = newtext + '\n'
        for key in self.info_label_strings_order:
          newtext = newtext + self.info_label_strings[key]
        
      self.info_label.text = newtext
    except:
      pass
  

  ##################################      
  ### SIMULATOR QUERYING METHODS ###
  ##################################
  def get_world_size(self):
    """ Returns the longest extent of the world/simulation box
    """    
    world_shape = self.sim.world.get_structure(self.sim.world.get_def_structure_id()).shape
    world_size  = 2.0 * world_shape.half_extent

    return max(world_size)
    
    
  def get_particle_data(self):
    
    particles, colors = [], []

    # Get particle data from simulator.
    for species in self.sim.world.species:
        for particle_id in self.sim.world.get_particle_ids(species):
          
            particle = self.sim.world.get_particle(particle_id)[1]
            particles.append(particle)
            
            try:
                color = self.particle_color_dict[species.id.serial]
            except:
                color = species.id.serial
                
            colors.append(color)

    return self.process_spheres(particles, colors)

        
  def get_protective_domain_data(self):
  
    spheres, sphere_colors = [], []
    cylinders, cylinder_colors = [], []

    if self.SHOW_DOMAINS == True:
      
        number_of_shells = self.sim.scheduler.size

        if number_of_shells > 0:
            top_event_id = self.sim.scheduler.top[0]
        else:
            top_event_id = None
         
        for object in self.sim.domains.itervalues():
          
            # First determine the type of shell 
            # and set the color accordingly
            
            # Single and pairs
            color = self.shell_color_dict[object.multiplicity]

            # Highlight top_event for singles and pairs
            if top_event_id != None and object.event_id == top_event_id:                
                color = self.shell_color_dict[0]                

            # Multi shells
            if isinstance(object, Multi):
                color = self.shell_color_dict[3]

            # Now get the shell data
            try:
                shell = object.shell_list[0][1]
                # Only cylinders have half_length.
                shell.shape.half_length                    
                cylinders.append(shell)
                cylinder_colors.append(color)                
            except:
                # Shell is not cylinder but sphere 
                # (Single, Pair or Multi)
                for _, shell in object.shell_list:
                    spheres.append(shell)
                    sphere_colors.append(color)

    return self.process_spheres(spheres, sphere_colors), \
              self.process_cylinders(cylinders, cylinder_colors)

              
  def get_cuboidal_region_data(self):
    
    boxes = [self.sim.world.get_structure(self.sim.world.get_def_structure_id()).shape]

    return self.process_boxes(boxes)

      
  def get_planar_surface_data(self):
    
    planes, colors = [], []
    
    world_id = self.sim.world.get_def_structure_id()
    planes = [surface for surface
                      in self.sim.world.structures
                      if isinstance(surface.shape, Plane)
                      and not surface.id == world_id]  # why this here? The world is not a plane!

    # Also make a list only containing the 'shape' objects (geometry) of the planes
    # We have to pass this to 'process_boxes' below because the geometrical features can only
    # be accessed via surface.shape.[..] (strangely enough, not so for the cylinders, see below)
    plane_shapes = [surface.shape for surface in planes]    
    
    for pl in planes:
      
      try:
        color = self.structure_color_dict[pl.sid.serial]
      except:
        color = pl.sid.serial
              
      colors.append(color)
      
    # In the visualizer, planes are treated as boxes with "zero" z-extension
    return self.process_boxes(plane_shapes, colors)
    
    
  def get_cylindrical_surface_data(self):
      
    cylinders, colors = [], []
    
    cylinders = [surface for surface
                         in self.sim.world.structures
                         if isinstance(surface.shape, Cylinder)
                         or isinstance(surface.shape, Disk)]
                         
    for cyl in cylinders:
      
        try:
          color = self.structure_color_dict[cyl.sid.serial]
        except:
          color = cyl.sid.serial
                
        colors.append(color)

    return self.process_cylinders(cylinders, colors)

    
  def get_spherical_surface_data(self):
    # TODO so far unused
    
    spheres = [surface for surface
                       in self.sim.world.structures
                       if isinstance(surface.shape, Sphere)]
    return self.process_spheres(spheres)


  ######################################      
  ### DATA FORMAT PROCESSING METHODS ###
  ######################################
  def process_spheres(self, spheres=[], color_list=[]):
    # Returns 3 lists:
    #   - positions, as (x,y,z) vectors
    #   - radii
    #   - colors
    #
    position_list, radius_list = [], []      

    for sphere in spheres:
        # The following is necessary to distinguish
        # shells from surfaces
        try:
            pos = sphere.position              
            radius = sphere.radius
        except:
            pos = sphere.shape.position
            radius = sphere.shape.radius              
                  
        formated_position = ( pos[0], pos[1], pos[2] )          
        
        position_list.append(formated_position)
        radius_list.append(radius)

    return (position_list, radius_list, color_list)
      

  def process_cylinders(self, cylinders=[], color_list=[]):
    # Return 5 lists:
    #   - positions
    #   - orientations
    #   - radii
    #   - lengths
    #   - colors
    #
    position_list, orientation_list, radius_list, length_list = [], [], [], []

    for cylinder in cylinders:
        # Pre-process disks which are put into the same
        # list as cylinders but lack the half_length property
        try:
            if isinstance(cylinder.shape, Disk):
                disk = cylinder
                cylinder = self.get_dummy_cylinder()
                cylinder.position    = disk.shape.position
                cylinder.radius      = disk.shape.radius
                cylinder.unit_z      = disk.shape.unit_z
                cylinder.half_length = 0.0
        except:
            pass

        # The following is necessary to distinguish
        # shells from surfaces
        try:
            position = cylinder.position
            radius = cylinder.radius
            orientation = cylinder.unit_z
            half_length = cylinder.half_length
        except:
            position = cylinder.shape.position
            radius = cylinder.shape.radius
            orientation = cylinder.shape.unit_z
            half_length = cylinder.shape.half_length
          
        # In eGFRD the position points to the center of the cylinder,
        # in Vpython to one of the ends, so we have to correct for this
        position = position - half_length * orientation

        position_list.append(position)
        orientation_list.append(orientation)
        radius_list.append(radius)
        length_list.append(2.0*half_length)

    return (position_list, orientation_list, radius_list, length_list, color_list)

    
  def process_boxes(self, boxes=[], color_list=[]):
    
    orientation = []
    position_list, orientation_list, size_list = [], [], []

    for box in boxes:
        
        position = self.to_vector(box.position)
      
        try:
            dz = box.half_extent[2]
        except IndexError:
            # Planes don't have half_extent[2].
            dz = ALMOST_ZERO
                            
        orientation_ux = self.to_vector(box.unit_x) # we take the x-axis to be "length"
        orientation_uy = self.to_vector(box.unit_y) # we also want the other vectors to be able 
        orientation_uz = self.to_vector(box.unit_z) # to fully construct the box when needed
        size = (2.0 * box.half_extent[0], 2.0 * box.half_extent[1], 2.0 * dz)
        
        position_list.append(position)
        orientation_list.append( [orientation_ux, orientation_uy, orientation_uz] )
        size_list.append(size)

    return (position_list, orientation_list, size_list, color_list)
    

  #############################      
  ### VISUALIZATION METHODS ###
  #############################
  def draw_cuboidal_region(self, position_list, orientation_list, size_list):

    if self.SHOW_SIM_BOX == True:
      
      # TODO Check that lists all have the same length
      for n in range(len(position_list)):
      
          orientation_ux, _, _ = orientation_list[n]
          
          if not(self.WIREFRAME_BOX):
            box = v.box(pos=position_list[n], axis=orientation_ux, size=size_list[n], color=(1, 1, 1), opacity=self.SIM_BOX_OPACITY)
          else:
            box = self.wireframe_box(pos=position_list[n], orientation=orientation_list[n], size=size_list[n], color=(0, 0, 1), radius=self.WIREFRAME_RADIUS)
            
          self.objects.append(box)
  
  
  def draw_planar_surfaces(self, position_list, orientation_list, size_list, color_list):
    
    if self.SHOW_STRUCTURES == True:
      
      # TODO Check that lists all have the same length
      for n in range(len(position_list)):
        
          orientation_ux, _, _ = orientation_list[n]
          
          if not(self.WIREFRAME_PLANES):
            # Standard planes are treated as boxes with inifinitesimal length in one of the three directions
            plane = v.box(pos=position_list[n], axis=orientation_ux, size=size_list[n], color=color_list[n], opacity=self.PLANAR_SURFACES_OPACITY)
          else:
            plane = box = self.wireframe_box(pos=position_list[n], orientation=orientation_list[n], size=size_list[n], color=color_list[n], radius=self.WIREFRAME_RADIUS)
            
          self.objects.append(plane)
  
  
  def draw_cylindrical_surfaces(self, position_list, orientation_list, radius_list, length_list, color_list):
      
    if self.SHOW_STRUCTURES == True:
      
      # TODO Check that lists all have the same length
      for n in range(len(position_list)):

          outer_radius = self.ROD_RADIUS_FACTOR*radius_list[n]
      
          if self.RODS_AS_HELICES == True and length_list[n] > 0.0:
            # The standard VPython helix object cannot be properly erased from the scene, so we use an own one ;-)            
            obj = self.tomeks_double_helix(pos=position_list[n], axis=orientation_list[n], radius=outer_radius, \
                                              length=length_list[n], color=color_list[n], thickness=self.HELIX_THICKNESS_FACTOR*outer_radius)
                                              # thickness=self.ROD_RADIUS_FACTOR*radius_list[n]/7 also worked well
          elif length_list[n] > 0.0:
            obj = v.cylinder(pos=position_list[n], axis=orientation_list[n], radius=outer_radius, \
                                  length=length_list[n], color=color_list[n], opacity=self.CYLINDER_SURFACES_OPACITY)
          else: # the cylinder actually is a disk
            obj = v.ring(pos=position_list[n], axis=orientation_list[n], radius=outer_radius, \
                                  thickness=self.HELIX_THICKNESS_FACTOR*outer_radius, color=color_list[n], opacity=self.DISK_SURFACES_OPACITY)
            
          self.objects.append(obj)


  def draw_particles(self, position_list, radius_list, color_list):  
    
    # TODO Check that lists all have the same length
    for n in range(len(position_list)):
      
        sphere = v.sphere(pos=position_list[n], radius=self.PARTICLE_RADIUS_FACTOR*radius_list[n], color=color_list[n])
        self.objects.append(sphere)
        
        
  def draw_protective_domains(self, spherical_shell_data, cylindrical_shell_data):

    if self.SHOW_DOMAINS == True:      
      
      # First unpack the data lists
      # TODO Check that lists all have the same length      
      Spos_list, Srad_list, Scol_list                       = spherical_shell_data
      Cpos_list, Cori_list, Crad_list, Clen_list, Ccol_list = cylindrical_shell_data
          
      # Draw the spherical shells
      for n in range(len(Spos_list)):
          
          used_radius  = Srad_list[n]
          used_opacity = self.SHELL_OPACITY
          
          # Special case: check whether the shell is a Multi (by checking its unique color)
          if Scol_list[n] == self.shell_color_dict[3]:
            # Then we want to scale up the radius by the same factor as the particles and make them less opaque
            used_radius  = self.PARTICLE_RADIUS_FACTOR * Srad_list[n]
            used_opacity = 0.5
          
          shell = v.sphere(pos=Spos_list[n], radius=used_radius, color=Scol_list[n], opacity=used_opacity)
          self.objects.append(shell)
          
      # Draw the cylindrical shells
      for n in range(len(Cpos_list)):
          
          # TODO can we somehow also scale up the radius of the cyl. shells of cyl. bound particles (exclusively!)
          shell = v.cylinder(pos=Cpos_list[n], axis=Cori_list[n], radius=Crad_list[n], length=Clen_list[n], color=Ccol_list[n], opacity=self.SHELL_OPACITY)
          self.objects.append(shell)

          
  def wireframe_box(self, pos=(0,0,0), orientation=[(1,0,0),(0,1,0),(0,0,1)], size=(1,1,1), color=(0,0,0), radius=0.1):
    
      # Make the input VPython vectors, such that we can use VPython's vector arithmetics
      P0 = v.vector(pos)
      # Normalize the orientatoin vectors
      Ux = v.norm(orientation[0])
      Uy = v.norm(orientation[1])
      Uz = v.norm(orientation[2])
      
      # Create a vector that holds the half-lengths of the box
      H = v.vector( (size[0]/2.0, size[1]/2.0, size[2]/2.0) )

      # Create the corner points
      # Nomenclature:
      #
      # - B = back side, F = front side
      # - LL = lower left, UL = upper left
      # - LR = lower right, UR = upper right
      #
      # All seen with respect to the 'axis' vector (Ux)
      #
      C_BLL = P0 - H[0]*Ux - H[1]*Uy - H[2]*Uz      
      C_BLR = P0 - H[0]*Ux + H[1]*Uy - H[2]*Uz
      C_BUL = P0 - H[0]*Ux - H[1]*Uy + H[2]*Uz
      C_BUR = P0 - H[0]*Ux + H[1]*Uy + H[2]*Uz
      C_FLL = P0 + H[0]*Ux - H[1]*Uy - H[2]*Uz      
      C_FLR = P0 + H[0]*Ux + H[1]*Uy - H[2]*Uz
      C_FUL = P0 + H[0]*Ux - H[1]*Uy + H[2]*Uz
      C_FUR = P0 + H[0]*Ux + H[1]*Uy + H[2]*Uz
      
      # Create lists that will hold the points for the separate curves
      curve_points_back  = [C_BLL, C_BUL, C_BUR, C_BLR, C_BLL]
      curve_points_front = [C_FLL, C_FUL, C_FUR, C_FLR, C_FLL]
      curve_points_left  = [C_BLL, C_BUL, C_FUL, C_FLL, C_BLL]
      curve_points_right = [C_BLR, C_BUR, C_FUR, C_FLR, C_BLR]
      
      # Now create a frame object and the separate frames, and associate the latter with the former
      f = v.frame()
      
      curve_back  = v.curve(pos=curve_points_back,  color=color, radius=radius, frame=f)
      curve_front = v.curve(pos=curve_points_front, color=color, radius=radius, frame=f)
      curve_left  = v.curve(pos=curve_points_left,  color=color, radius=radius, frame=f)
      curve_right = v.curve(pos=curve_points_right, color=color, radius=radius, frame=f)
      
      return f
  
  
  def wireframe_plane(self, pos=(0,0,0), orientation=[(1,0,0),(0,1,0),(0,0,1)], size=(1,1,1), color=(0,0,0), radius=0.1):
    
      # Make the input VPython vectors, such that we can use VPython's vector arithmetics
      P0 = v.vector(pos)
      # Normalize the orientatoin vectors
      Ux = v.norm(orientation[0])
      Uy = v.norm(orientation[1])
      #Uz = v.norm(orientation[2]) # not used here
      
      # Create a vector that holds the half-lengths of the box
      H = v.vector( (size[0]/2.0, size[1]/2.0, 0.0) )

      # Create the corner points
      # Nomenclature:
      #
      # - BL = back left, BR = back right
      # - FL = front left, FR = front right
      #
      # All seen with respect to the 'axis' vector (Ux)
      C_BL = P0 - H[0]*Ux - H[1]*Uy
      C_BR = P0 - H[0]*Ux + H[1]*Uy      
      C_FL = P0 + H[0]*Ux - H[1]*Uy
      C_FR = P0 + H[0]*Ux + H[1]*Uy
      
      # Create lists that will hold the points for the curve
      curve_points  = [C_BL, C_FL, C_FR, C_BR, C_BL]
      
      # Now create the curve object      
      curva_yebana  = v.curve(pos=curve_points, color=color, radius=radius)
      
      return curva_yebana

      
  def tomeks_helix(self, pos=(0,0,0), axis=(1,0,0), radius=1.0, length=10.0, color=(1,0,0), thickness=0.01):
          
      Ncoils  = self.N_HELIX_COILS
      Npoints = 10 * Ncoils 
      
      # Draw the helix in the default coordinate system and tie it to a frame
      f = v.frame()
      c = v.curve( x = v.arange(0,length,length/Npoints), frame=f, color=color, radius = thickness)
      c.y = radius * np.sin(c.x * 2.0*np.pi/length * Ncoils)
      c.z = radius * np.cos(c.x * 2.0*np.pi/length * Ncoils)
      
      # Now rotate it to the right orientation
      f.pos  = pos
      f.axis = axis
      
      return f
      
      
  def tomeks_double_helix(self, pos=(0,0,0), axis=(1,0,0), radius=1.0, length=10.0, color=(1,0,0), thickness=0.01):

      # This is just the two of the standard 'tomeks_helices' with an offset of 1/2 of a coil
      Ncoils  = self.N_HELIX_COILS
      Npoints = 10 * Ncoils 
      
      offset = length / Ncoils / 2.0
      
      # Draw the helices in the default coordinate system and tie it to a frame
      f = v.frame()
      
      c1 = v.curve( x = v.arange(0,length,length/Npoints), frame=f, color=color, radius = thickness)
      c1.y = radius * np.sin(c1.x * 2.0*np.pi/length * Ncoils)
      c1.z = radius * np.cos(c1.x * 2.0*np.pi/length * Ncoils)
      
      c2 = v.curve( x = v.arange(0,length,length/Npoints), frame=f, color=self.complementary(color), radius = thickness)
      c2.y = radius * np.sin( (c2.x - offset) * 2.0*np.pi/length * Ncoils)
      c2.z = radius * np.cos( (c2.x - offset) * 2.0*np.pi/length * Ncoils)
      
      # Now rotate it to the right orientation
      f.pos  = pos
      f.axis = axis
      
      return f
      

  ######################      
  ### HELPER METHODS ###
  ######################
  def complementary(self, RGB_color):
    
      # Get the/a complementary color of a RGB color (passed as a tuple)      
      # Convert to HSV system first; then we can change the hue
      HSV_color = colorsys.rgb_to_hsv(RGB_color[0], RGB_color[1], RGB_color[2])            
      H = HSV_color[0]      
      S = HSV_color[1]
      V = HSV_color[2]
      
      # Invert the hue (turn the 'hue wheel' by 120 degrees)
      H = (H + 0.5) % 1.0
      
      return colorsys.hsv_to_rgb(H,S,V)
  
  def to_vector(self, v):
  
    # TODO Assert v has only 3 entries    
    return ( v[0], v[1], v[2] )

    
  def normalize_vector(self, v):
    
    norm = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])    
    nv   = (v[0]/norm, v[1]/norm, v[2]/norm)
    
    return nv
    
    
  def get_dummy_sphere(self):
      return Sphere([0, 0, 0], ALMOST_ZERO)

      
  def get_dummy_cylinder(self):
      return Cylinder([0, 0, 0], ALMOST_ZERO, [0, 0, 1], ALMOST_ZERO)

      
  def get_dummy_box(self):
      return Box([0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1],
                  ALMOST_ZERO, ALMOST_ZERO, ALMOST_ZERO)

  def cisleithanian_name(self):
      return "Der k. u. k. eGFRD-Visualisator"
      
  def transleithanian_name(self):
      return "Az csaszari es kiralyi eGFRD vizualizacioja szoftver"



# A derived / wrapper class with a shorter (i.e. usable) name
class TIRViE(TheImperialRoyalVisualizerForEGFRD):
  """

    This is just an alias for 'TheImperialRoyalVisualizerForEGFRD'
  
  """
  def __init__(self, eGFRDsim, particle_color_dict={}, structure_color_dict={}, capture=False, capture_dir='./tirvie_frames'):

    TheImperialRoyalVisualizerForEGFRD.__init__(self,eGFRDsim, particle_color_dict, structure_color_dict, capture, capture_dir)
    
    