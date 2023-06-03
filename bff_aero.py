import h5py as h5
import numpy as np


class BFF_Aero:
    """
        BFF_Aero contains all attributes to define the aerodynamic grid and makes it accessible for SHARPy.
    """
    def __init__(self, m, structure, case_name, case_route, **kwargs):
        """        

        Key-Word Arguments:
            - m: number of chordwise vortex ring panels
            - structure: structural object of BFF model
        """
        self.m = m
        self.structure = structure

        self.route = case_route
        self.case_name = case_name

        self.elastic_axis_wing = 0.5
        self.elastic_axis_mid_body = abs(self.structure.x_nose / self.structure.chord_mid_body)
        

    def generate(self):
        """
            Function to set up all necessary parameter inputs to define the geometry and discretisation
            of the lifting surfaces. Finally, informations are saved to an .h5 file serving as an input
            file for SHARPy.
        """
        self.set_wing_properties()
        self.set_control_surfaces()
        self.mirror_wing()
        self.write_input_file()

    def set_wing_properties(self):
        """
            Sets necessary parameters to define the lifting surfaces of one wing (right).
        """
        self.airfoil_distribution = np.zeros((self.structure.n_elem, self.structure.n_node_elem), dtype=int)
        self.surface_distribution = np.zeros((self.structure.n_elem,), dtype=int)
        self.surface_m = np.zeros((self.structure.n_surfaces, ), dtype=int) + self.m

        self.aero_node = np.zeros((self.structure.n_node,), dtype=bool)
        self.aero_node[:] = True

        self.twist = np.zeros((self.structure.n_elem, self.structure.n_node_elem))
        self.sweep = np.zeros_like(self.twist)
        self.chord = np.zeros_like(self.twist)
        self.elastic_axis = np.zeros_like(self.twist)

        #wing
        self.chord[self.structure.n_elem_body:self.structure.n_elem//2] = self.structure.chord_wing
        self.elastic_axis[self.structure.n_elem_body:self.structure.n_elem//2] = self.elastic_axis_wing 
        # body
        for i_elem in range(self.structure.n_elem_body):
            for i_local_node in range(self.structure.n_node_elem):
                inode = self.structure.conn[i_elem, i_local_node]
                self.chord[i_elem, i_local_node], self.elastic_axis[i_elem, i_local_node] = self.get_chord_and_elastic_axis_body_node(inode)

    def mirror_surface(self):
        """
            Mirrors the parameters from the right lifting surface for the left one.
        """
        self.elastic_axis[self.structure.n_elem_right:, :] = self.elastic_axis[:self.structure.n_elem_right, :]
        self.chord[self.structure.n_elem_right:, :] = self.chord[:self.structure.n_elem_right, :]
        self.surface_distribution[self.structure.n_elem_right:] += 1 
    
    def get_chord_and_elastic_axis_body_node(self, inode):
        """
            Calculates the chord and elastic of the lifting surface of the body at each structural node 
            based on defined geometry planform.
        """  
        delta_x_TE = np.tan(self.structure.sweep_body_TE) * self.structure.y[inode]
        if self.structure.y[inode] <= self.structure.halfspan_body / 2.:
            delta_x_LE =  np.tan(self.structure.sweep_body_LE) * self.structure.y[inode]
        else:
            delta_x_LE = self.structure.offset_nose_quarter_body_LE + np.tan(self.structure.sweep_quarter_body_LE) * (self.structure.y[inode] - self.structure.halfspan_body / 2.)
        chord =  self.structure.chord_mid_body - delta_x_LE - delta_x_TE
        elastic_axis = (abs(self.structure.x_nose) - delta_x_LE) / chord
        return chord, elastic_axis

    
    def write_input_file(self):
        """
            Writes previously defined parameters to an .h5 file which serves later as an 
            input file for SHARPy.
        """
            
        with h5.File(self.route + '/' + self.case_name + '.aero.h5', 'a') as h5file:
            airfoils_group = h5file.create_group('airfoils')
            # add one airfoil
            airfoils_group.create_dataset('0', 
                                          data=np.column_stack(
                                            self.generate_naca_camber(P=0, M=0)
                                          ))

            h5file.create_dataset('chord', data=self.chord)
            h5file.create_dataset('twist', data=self.twist)
            h5file.create_dataset('sweep', data=self.sweep)

            # airfoil distribution
            h5file.create_dataset('airfoil_distribution', data=self.airfoil_distribution)
            h5file.create_dataset('surface_distribution', data=self.surface_distribution)
            h5file.create_dataset('surface_m', data=self.surface_m)
            h5file.create_dataset('m_distribution', data='uniform')
            h5file.create_dataset('aero_node', data=self.aero_node)
            h5file.create_dataset('elastic_axis', data=self.elastic_axis)
            h5file.create_dataset('control_surface', data=self.control_surface)
            h5file.create_dataset('control_surface_deflection', data=self.control_surface_deflection)
            h5file.create_dataset('control_surface_chord', data=self.control_surface_chord)
            h5file.create_dataset('control_surface_hinge_coord', data=self.control_surface_hinge_coord)
            h5file.create_dataset('control_surface_type', data=self.control_surface_type)
   

    def generate_naca_camber(self,M=0, P=0):
        """
            Generates the camber line coordinates of a specified NACA airfoil.
        """
        mm = M*1e-2
        p = P*1e-1

        def naca(x, mm, p):
            if x < 1e-6:
                return 0.0
            elif x < p:
                return mm/(p*p)*(2*p*x - x*x)
            elif x > p and x < 1+1e-6:
                return mm/((1-p)*(1-p))*(1 - 2*p + 2*p*x - x*x)

        x_vec = np.linspace(0, 1, 1000)
        y_vec = np.array([naca(x, mm, p) for x in x_vec])
        return x_vec, y_vec
