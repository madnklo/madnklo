import os
import shutil
root_path = os.path.split(os.path.dirname(os.path.realpath( __file__ )))[0]
pjoin = os.path.join

class CurrentsExporterException(Exception):
    pass

def copy_while_creating_directories(srcpath, dstpath):
    base_directory = os.path.split(dstpath)[0]
    if not os.path.isdir(base_directory):
        os.makedirs(base_directory)
    if os.path.isdir(srcpath):
        shutil.copytree(srcpath,dstpath)
    else:
        shutil.copy(srcpath,dstpath)

class GenericCurrentsExporter(object):
    """ Simple current exporter that only takes care of copying current source files
    to the intended destination."""

    def __init__(self, relative_resource_paths=None):
        """ Stores the list of ressource files/directories to copy to the destination for that
        particular subtraction current scheme."""

        # If not specified then copy everything
        if relative_resource_paths is None: 
            self.relative_resource_paths = \
                [ path for path in os.listdir(root_path) if 
                  (not any(path.endswith(veto_extension)) for veto_extension in ['.pyc']) ]
        else:
            self.relative_resource_paths = relative_resource_paths
            mandatory_resources = ['__init__.py', 'commons', 'subtraction_schemes/__init__.py']
            for mandatory_resource in mandatory_resources:
                if mandatory_resource not in self.relative_resource_paths:
                    self.relative_resource_paths.append(mandatory_resource)

        # Now specify rust-specific paths
        self.rust_template_files = {
            'all_subtraction_current_evaluators.rs' : pjoin('commons','rust','all_subtraction_current_evaluators.rs'),
            'subtraction_current_evaluator.rs': pjoin('commons', 'rust', 'subtraction_current_evaluator.rs'),
            'subtraction_current.rs': pjoin('commons', 'rust', 'subtraction_current.rs'),
        }

    def export_rust(self, output_directory):
        """ Exports rust *common* resources. The ones for each subtraction current will be implemented via
        dynamic meta-coding."""

        for file_name in ['subtraction_current.rs']:
            shutil.copy(self.rust_template_files[file_name],pjoin(output_directory,file_name))


    def export(self, destination_root_path, mapped_currents_to_export=None, **opts):
        """ Copy all resources to their destination path. """
        
        # For this exporter there is no differentiation as to what needs to be exported
        # and we simply export all relative_resource_paths specified.
        # This therefore needs to be done once only irrespectively of which currents
        # needs to be exported.
        if os.path.isdir(destination_root_path):
            return

        # First make sure that the target directory to copy things into exists
        output_base_path = os.path.split(destination_root_path)[0]
        if not os.path.isdir(output_base_path):
            raise CurrentsExporterException( "The specified destination path "+
                "'%s' for exporting subtraction current cannot be found."%output_base_path )

        for path in self.relative_resource_paths:
            copy_while_creating_directories(pjoin(root_path, path),
                                            pjoin(destination_root_path, path) )
