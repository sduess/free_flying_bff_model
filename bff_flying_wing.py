import configobj
from bff_structure import BFF_Structure
from bff_aero import BFF_Aero
import os
import sharpy.sharpy_main


class BFF_Flying_Wing:

    def __init__(self, case_name, case_route, output_route):
        self.case_name = case_name
        self.case_route = case_route
        self.output_route = output_route

        self.structure = None
        self.aero = None

        self.settings = None

    def init_aeroelastic(self, **kwargs):
        self.clean()
        self.init_structure(**kwargs)
        self.init_aero(**kwargs)

    def init_structure(self, **kwargs):
        self.structure = BFF_Structure(self.case_name, self.case_route, **kwargs)

    def init_aero(self, m, **kwargs):
        self.aero = BFF_Aero(m, self.structure, self.case_name, self.case_route,**kwargs)

    def generate(self):

        if not os.path.isdir(self.case_route):
            os.makedirs(self.case_route)

        self.structure.generate()

        if self.aero is not None:
            self.aero.generate()


    def create_settings(self, settings):
        file_name = self.case_route + '/' + self.case_name + '.sharpy'
        config = configobj.ConfigObj()
        config.filename = file_name
        for k, v in settings.items():
            config[k] = v
        config.write()
        self.settings = settings

    def clean(self):
        list_files = ['.fem.h5', '.aero.h5', '.nonlifting_body.h5', '.dyn.h5', '.mb.h5', '.sharpy', '.flightcon.txt']
        for file in list_files:
            path_file = self.case_route + '/' + self.case_name + file
            if os.path.isfile(path_file):
                os.remove(path_file)

    def run(self):
        sharpy.sharpy_main.main(['', self.case_route + '/' + self.case_name + '.sharpy'])

