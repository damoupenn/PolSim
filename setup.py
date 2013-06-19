import glob

__version__='1.0.0'

setup_args = {
        'name': 'PolSim',
        'author': 'David F. Moore',
        'author_email': 'damo@sas.upenn.edu',
        'license': 'GPL',
        'packages': ['PolSim'],
        'package_dir': {'PolSim':'src'},
        'package_data': {'data':['input_data/*'],
                            'config':['config_files/*']},
        'scripts': glob.glob('scripts/*.*'),
        'version': __version__,
        }

if __name__ == "__main__":
    from distutils.core import setup
    apply(setup, (), setup_args)
