import h5py as h5
import numpy as np


class BFF_Structure:

    """
        BFF_Structure contains all attributes to define the beam geometriy 
        and discretisation and makes it accessible for SHARPy.
    """
    def __init__(self, case_name, case_route, **kwargs):
        self.sigma = kwargs.get('sigma', 1)
        self.n_elem_multiplier = kwargs.get('n_elem_multiplier', 1)

        self.route = case_route
        self.case_name = case_name

        self.n_elem = None
        self.n_node = None
        self.n_node_elem = 3
        self.n_surfaces = 2

        self.x = None
        self.y = None
        self.z = None

        self.n_elem_body = None
        self.n_elem_wing = None
    
        self.n_node_body = None
        self.n_node_wing = None

        # TODO: store as aircraft information?
        self.chord_wing = 0.0957 # m
        self.chord_mid_body = 0.2743 # m
        self.chord_quarter_body = 0.1446 # m

        self.halfspan_wing = 0.3786 + 0.0607 * 2 # m
        self.halfspan_body = 0.0607 * 2 #m

        self.offset_nose_quarter_body_LE = 0.1052
        self.offset_nose_wing_start_LE = 0.1295      
        self.x_offset_nose_tip_LE = 0.2673 
        self.x_nose = -(self.chord_wing / 2. + self.offset_nose_wing_start_LE)
        self.set_sweep_angles()

        self.thrust = 0.

    def generate(self):
        """
            Function to set up all necessary parameter inputs to define the geometry and discretisation
            of the beam. Finally, informations are saved to an .h5 file serving as an input
            file for SHARPy.
        """
        self.set_element_and_nodes()

        self.set_stiffness_and_mass_propoerties()
        self.set_beam_properties()


        self.set_lumped_masses()
        self.mirror_beam()

        self.app_forces = np.zeros((self.n_node, 6))        
        self.app_forces[0] = [0, self.thrust, 0, 0, 0, 0]
        self.write_input_file()



    def set_element_and_nodes(self):
        self.n_elem_body =  2*int(2*self.n_elem_multiplier)
        self.n_elem_wing =  3*int(1 * self.n_elem_body)
        self.n_elem = 2 * (self.n_elem_body + self.n_elem_wing)

        self.n_node_body = self.n_elem_body*(self.n_node_elem - 1) + 1
        self.n_node_wing = self.n_elem_wing*(self.n_node_elem - 1)
        self.n_node = 2 * (self.n_node_body + self.n_node_wing) - 1

        self.n_node_right = self.n_node_body + self.n_node_wing
        self.n_elem_right = self.n_elem // 2

    def set_sweep_angles(self):
        """
            Calculates all sweep angles of leading and trailing edges for later use.
        """
        # wing beam
        delta_y_tip = self.halfspan_wing - self.halfspan_body
        delta_x_tip = self.x_offset_nose_tip_LE - self.offset_nose_wing_start_LE
        self.sweep_wing = np.arctan(delta_x_tip / delta_y_tip)

        # body LE
        delta_x_mid_body_LE = self.offset_nose_quarter_body_LE 
        self.sweep_body_LE = np.arctan(delta_x_mid_body_LE / (self.halfspan_body/2))

        # body TE 
        delta_x_mid_body_TE = self.chord_mid_body - self.offset_nose_wing_start_LE - self.chord_wing
        self.sweep_body_TE = np.arctan(delta_x_mid_body_TE / self.halfspan_body)

        # quarter body LE 
        delta_x_quarter_body_LE = self.offset_nose_wing_start_LE - self.offset_nose_quarter_body_LE 
        self.sweep_quarter_body_LE = np.arctan(delta_x_quarter_body_LE / (self.halfspan_body/2))

    def set_lumped_masses(self):
        # TODO: add accelerometer lumped masses
        n_lumped_mass = 1
        self.lumped_mass_nodes = np.zeros((n_lumped_mass, ), dtype=int)
        self.lumped_mass = np.zeros((n_lumped_mass, ))
        self.lumped_mass_inertia = np.zeros((n_lumped_mass, 3, 3))
        self.lumped_mass_position = np.zeros((n_lumped_mass, 3))
        # lumped mass to get mass properties of wing
        self.lumped_mass[0] = 0.16 - 0.07235312427332974
        self.lumped_mass_position[0, 1] = 0.053474066

    def set_beam_coordinate_nodes(self):        
        """
            Defines the coordinates of each beam node.
        """     
        self.x = np.zeros((self.n_node, ))
        self.y = np.zeros((self.n_node, ))
        self.z = np.zeros((self.n_node, ))
        
        self.y[:self.n_node_body] = np.linspace(0.0, self.halfspan_body, self.n_node_body)
        self.y[self.n_node_body:self.n_node_right] = np.linspace(self.halfspan_body, self.halfspan_wing, self.n_node_wing + 1)[1:]
        # account for beam sweep in x coordinate
        self.x[self.n_node_body:self.n_node_body + self.n_node_wing] = (self.y[self.n_node_body:self.n_node_body + self.n_node_wing]-self.halfspan_body) * np.tan(self.sweep_wing)
        
    def set_stiffness_and_mass_propoerties(self):
        
        n_material = 3 # body + wing right + wing left
        self.stiffness = np.zeros((n_material, 6, 6))
        self.mass = np.zeros((n_material, 6, 6))
        self.elem_stiffness = np.zeros((self.n_elem, ), dtype=int)
        self.mass = np.zeros((n_material, 6, 6))
        self.elem_mass = np.zeros((self.n_elem, ), dtype=int)
        
        # TODO: add inertia to (mass???)
        ea = 1e7
        ga = 1e5 
        gj = 1e4
        eiy = 1e4 
        eiz = eiy
        # gj = 0.06 # Torsional stiffness
        # ga = gj * 20 # Shear stiffness in the local y axis
        # eiy = 0.18 # Bending stiffness around the flapwise direction
        # eiz = 0.18 # Bending stiffness around the edgewise direction
        # ea = eiy * 1e4 # Axial stiffness
        m_bar_main = 0.069 
        j_bar_main = m_bar_main / 10.
        base_stiffness = self.sigma * np.diag([ea, ga, ga, gj, eiy, eiz]) / 1e4
        base_mass = np.diag([m_bar_main, m_bar_main, m_bar_main, j_bar_main, 0.5 * j_bar_main, 0.5 * j_bar_main])

        sigma_main_body = 100

        # Stiffness and mass properties
        self.stiffness[0, ...] = base_stiffness * sigma_main_body
        self.stiffness[1, ...] = base_stiffness
        self.stiffness[2, ...] = base_stiffness

        self.mass[:, ...] =  base_mass

    def set_beam_properties(self):
        """
            Defines all necessary parameters to define the beam including node coordinate,
            elements with their associated nodes, frame of reference delta, boundary conditions
            such as reference node and free tips, stiffness and mass property ID for each element,
            and twist.
        """     
        self.set_beam_coordinate_nodes()
        
        self.frame_of_reference_delta = np.zeros((self.n_elem, self.n_node_elem, 3))
        self.conn = np.zeros((self.n_elem, self.n_node_elem), dtype=int)
        for ielem in range(self.n_elem // 2):
            self.conn[ielem, :] = ((np.ones((3, )) * ielem * (self.n_node_elem - 1)) +
                                [0, 2, 1])               
            for ilocalnode in range(self.n_node_elem):
                self.frame_of_reference_delta[ielem, ilocalnode, :] = [-1.0, 0.0, 0.0]  
            if ielem > self.n_elem_body:
                # wing and body properties differ
                self.elem_stiffness[ielem] += 1
                self.elem_mass[ielem] += 1
        
        self.beam_number = np.zeros((self.n_elem, ), dtype=int)
        self.boundary_conditions = np.zeros((self.n_node, ), dtype=int)
        self.boundary_conditions[0] = 1
        self.boundary_conditions[self.n_node_wing + self.n_node_body] = -1 # free tip

        self.structural_twist = np.zeros((self.n_elem, self.n_node_elem))

    def calculate_aircraft_mass(self):
        # get structural mass for each component (beam ID)
        list_elem_mass = []
        center_of_gravity = np.zeros((3, ))
        for i_elem in range(self.n_elem):
            start_node = self.conn[i_elem, 0]
            end_node = self.conn[i_elem, 1]
            # calculate length assuming that elem is straight (unloaded)
            length_elem = np.sqrt((self.x[start_node]-self.x[end_node])**2 
                                  + (self.y[start_node]-self.y[end_node])**2 
                                  +(self.z[start_node]-self.z[end_node])**2 )
            distance = [self.x[self.conn[i_elem, -1]],
                        self.y[self.conn[i_elem, -1]],
                        self.z[self.conn[i_elem, -1]]]
            distributed_mass_elem = self.mass[self.elem_mass[i_elem], 0, 0]
            mass_elem = distributed_mass_elem * length_elem
            list_elem_mass.append(mass_elem)
            
            for i_dim in range(3):
                center_of_gravity[i_dim] += mass_elem * distance[i_dim]
        total_mass_structure = sum(list_elem_mass)
        
        print("Total structural mass = ", total_mass_structure)
        for i_beam in set(self.beam_number):
            structural_mass_beam = sum(np.array(list_elem_mass)[self.beam_number == int(i_beam)])
            print("Total structural mass for beam {} is {} kg".format(i_beam, structural_mass_beam))
        
        # Get lumped masses
        n_lumped_masses_wing = len(self.lumped_mass)
        for i_mass in range(n_lumped_masses_wing):      
            position_B_frame = self.lumped_mass_position[i_mass, :]
            position_G_frame = np.zeros_like(position_B_frame)
            position_G_frame[0] = -position_B_frame[1]
            position_G_frame[1] = position_B_frame[0]
            position_G_frame[2] = position_B_frame[2]
                
            if self.lumped_mass[i_mass] > 0:
                center_of_gravity[0] +=  self.lumped_mass[i_mass] * (position_G_frame[0] + self.x[self.lumped_mass_nodes[i_mass]])
                center_of_gravity[1] +=  self.lumped_mass[i_mass] * (position_G_frame[1] + self.y[self.lumped_mass_nodes[i_mass]])
                center_of_gravity[2] +=  self.lumped_mass[i_mass] * (position_G_frame[2] + self.z[self.lumped_mass_nodes[i_mass]])

        total_mass_lumped_masses = sum(self.lumped_mass)
        total_mass = total_mass_lumped_masses + total_mass_structure
        center_of_gravity /= total_mass
   
        print("Total lumped masses ", total_mass_lumped_masses)
        print("Total mass aircraft = ", total_mass)
        print("Center of Gravity = ", center_of_gravity)

    def mirror_beam(self):
        """
            Mirrors the parameters from the beam representing the right free-flying wing
            for the left one.
        """
        self.x[self.n_node_right:]  = self.x[1:self.n_node_right]
        self.y[self.n_node_right:]  = - self.y[1:self.n_node_right]
        self.frame_of_reference_delta[self.n_elem_right:, :, :] = self.frame_of_reference_delta[:self.n_elem_right, :, :] * (-1)
        self.elem_stiffness[self.n_elem_right:] = self.elem_stiffness[:self.n_elem_right] 
        self.elem_mass[self.n_elem_right:] = self.elem_mass[:self.n_elem_right]
        # wing might have different properties eventually
        self.elem_stiffness[self.n_elem_right + self.n_elem_body + 1:] = 2
        self.elem_mass[self.n_elem_right + self.n_elem_body + 1:] = 2
        self.beam_number[self.n_elem_right:] = 1        
        self.boundary_conditions[-1] = -1 # free tip
        self.conn[self.n_elem_right:, :] = self.conn[:self.n_elem_right, :] + self.n_node_right - 1
        self.conn[self.n_elem_right, 0] = 0 

    def set_thrust(self, value):
        self.thrust = value

    def write_input_file(self):
        """
            Writes previously defined parameters to an .h5 file which serves later as an 
            input file for SHARPy.
        """                
        with h5.File(self.route + '/' + self.case_name + '.fem.h5', 'a') as h5file:
            h5file.create_dataset('coordinates', data=np.column_stack((self.x, self.y, self.z)))
            h5file.create_dataset('connectivities', data=self.conn)
            h5file.create_dataset('num_node_elem', data=self.n_node_elem)
            h5file.create_dataset('num_node', data=self.n_node)
            h5file.create_dataset('num_elem', data=self.n_elem)
            h5file.create_dataset('stiffness_db', data=self.stiffness)
            h5file.create_dataset('elem_stiffness', data=self.elem_stiffness)
            h5file.create_dataset('mass_db', data=self.mass)
            h5file.create_dataset('elem_mass', data=self.elem_mass)
            h5file.create_dataset('frame_of_reference_delta', data=self.frame_of_reference_delta)
            h5file.create_dataset('boundary_conditions', data=self.boundary_conditions)
            h5file.create_dataset('beam_number', data=self.beam_number)

            h5file.create_dataset('structural_twist', data=self.structural_twist)
            h5file.create_dataset('app_forces', data=self.app_forces)
            h5file.create_dataset('lumped_mass_nodes', data=self.lumped_mass_nodes)
            h5file.create_dataset('lumped_mass', data=self.lumped_mass)
            h5file.create_dataset('lumped_mass_inertia', data=self.lumped_mass_inertia)
            h5file.create_dataset('lumped_mass_position', data=self.lumped_mass_position)