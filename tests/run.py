# ________________________________________________________________________________
#
# Script run.py
#
# Runs tests for CI
#
# ________________________________________________________________________________

# ________________________________________________________________________________
# Imports

import os
import importlib
import sys
import argparse
import subprocess
import json
import time

# ________________________________________________________________________________
# Parameters

# Configurations for running the tests
configuration_list = {

    'sequential' : { 'compiler' : None,
                 'cmake' : '-DCMAKE_VERBOSE_MAKEFILE=ON -DBACKEND="sequential"',
                 'env' : '',
                 'prefix' : '',
                 'args' : [ '', '', '', '', ''],
                 'exe_name' : 'minipic.seq',
                 'threads' : [8, 8, 8, 1, 1],
                 'benchmarks' : [ "default", "beam", "antenna", "Ecst", "Bcst" ]},

    'openmp' : { 'compiler' : None,
                 'cmake' : '-DCMAKE_VERBOSE_MAKEFILE=ON -DBACKEND="openmp"',
                 'env' : 'OMP_PROC_BIND=spread',
                 'prefix' : '',
                 'args' : [ '', '', '', '', ''],
                 'exe_name' : 'minipic.omp',
                 'threads' : [8, 8, 8, 1, 1],
                 'benchmarks' : [ "default", "beam", "antenna", "Ecst", "Bcst" ]},

    'kokkos' : { 'compiler' : 'clang++',
                 'cmake' : '-DCMAKE_VERBOSE_MAKEFILE=ON -DBACKEND="kokkos"',
                 'env' : 'OMP_PROC_BIND=spread',
                 'prefix' : '',
                 'args' : [ '', '', '', '-it 10000', '-it 10000'],
                 'exe_name' : 'minipic.kokkos',
                 'threads' : [8, 8, 8, 1, 1],
                 'benchmarks' : [ "default", "beam", "antenna", "Ecst", "Bcst" ]},

    'kokkos_gpu' : { 'compiler' : 'g++',
                 'cmake' : '-DCMAKE_VERBOSE_MAKEFILE=ON -DBACKEND="kokkos" -DDEVICE="nvidia_a100" -DCMAKE_BUILD_TYPE=Release',
                #  'cmake' : '-DCMAKE_VERBOSE_MAKEFILE=ON -DBACKEND="kokkos" -DDEVICE="nvidia_v100"',
                 'env' : 'OMP_PROC_BIND=spread',
                 'prefix' : '',
                 'args' : [ '', '', '', '', ''],
                 'exe_name' : 'minipic.kokkos',
                 'threads' : [1, 1, 1],
                 'benchmarks' : [ "default_gpu", "beam_gpu", "antenna"]},

    'kokkos_unified' : { 'compiler' : 'clang++',
                 'cmake' : '-DCMAKE_VERBOSE_MAKEFILE=ON -DBACKEND="kokkos_unified"',
                 'env' : 'OMP_PROC_BIND=spread',
                 'prefix' : '',
                 'args' : [ '', '', '', '-it 10000', '-it 10000'],
                 'exe_name' : 'minipic.kokkos_unified',
                 'threads' : [8, 8, 8, 1, 1],
                 'benchmarks' : [ "default", "beam", "antenna", "Ecst", "Bcst" ]},

    'kokkos_unified_gpu' : { 'compiler' : None,
                 #'cmake' : '-DCMAKE_VERBOSE_MAKEFILE=ON -DBACKEND="kokkos_unified" -DDEVICE="nvidia_v100"',
                 'cmake' : '-DCMAKE_VERBOSE_MAKEFILE=ON -DBACKEND="kokkos_unified" -DDEVICE="nvidia_a100"',
                 'env' : 'OMP_PROC_BIND=spread',
                 'prefix' : '',
                 'args' : [ '', '', '', '', ''],
                 'exe_name' : 'minipic.kokkos_unified',
                 'threads' : [1, 1, 1],
                 'benchmarks' : [ "default_gpu", "beam_gpu", "antenna"]},

    'kokkos_dualview_unified_gpu' : { 'compiler' : None,
                 #'cmake' : '-DCMAKE_VERBOSE_MAKEFILE=ON -DBACKEND="kokkos_dualview_unified" -DDEVICE="nvidia_v100"',
                 'cmake' : '-DCMAKE_VERBOSE_MAKEFILE=ON -DBACKEND="kokkos_dualview_unified" -DDEVICE="nvidia_a100"',
                 'env' : 'OMP_PROC_BIND=spread',   
                 'prefix' : '',
                 'args' : [ '', '', '', '', ''],
                 'exe_name' : 'minipic.kokkos_dualview_unified',
                 'threads' : [1, 1, 1],
                 'benchmarks' : [ "default_gpu", "beam_gpu", "antenna"]},

    'thrust_gpu' : { 'compiler' : 'nvcc',
                 'cmake' : '-DCMAKE_VERBOSE_MAKEFILE=ON -DBACKEND="thrust" -DDEVICE="nvidia_a100" -DCMAKE_BUILD_TYPE=Release',
                 'env' : '',
                 'prefix' : '',
                 'args' : [ '', '', '', '', ''],
                 'exe_name' : 'minipic.thrust',
                 'threads' : [1, 1, 1],
                 'benchmarks' : [ "default_gpu", "beam_gpu", "antenna"]},

    'openmp_task' : { 'compiler' : 'clang++',
                 'cmake' : '-DCMAKE_VERBOSE_MAKEFILE=ON -DBACKEND="openmp_task"',
                 'env' : 'OMP_PROC_BIND=spread OMP_MAX_ACTIVE_LEVELS=10',
                 'prefix' : '',
                 'args' : [ '', '', '', '', ''],
                 'exe_name' : 'minipic.omp_task',
                 'threads' : [8, 8, 8, 1, 1],
                 'benchmarks' : [ "default", "beam", "antenna", "Ecst", "Bcst" ]},

    'openmp_target' : { 'compiler' : 'nvc++',
                 'cmake' : '-DCMAKE_VERBOSE_MAKEFILE=ON -DBACKEND="openmp_target"',
                 'env' : 'OMP_PROC_BIND=spread OMP_MAX_ACTIVE_LEVELS=10',
                 'prefix' : '',
                 'args' : [ '', '', '', '', ''],
                 'exe_name' : 'minipic.omp_target',
                 'threads' : [1, 1, 1, 1, 1],
                 'benchmarks' : [ "default_gpu", "beam_gpu", "antenna" ]},

    'openacc' : { 'compiler' : 'nvc++',
                 #'cmake' : '-DCMAKE_VERBOSE_MAKEFILE=ON -DBACKEND="openacc" -DDEVICE="nvidia_v100"',
                 'cmake' : '-DCMAKE_VERBOSE_MAKEFILE=ON -DBACKEND="openacc" -DDEVICE="nvidia_a100"',
                 'env' : '',
                 'prefix' : '',
                 'args' : [ '', '', '', '', ''],
                 'exe_name' : 'minipic.acc',
                 'threads' : [1, 1, 1, 1, 1],
                 'benchmarks' : [ "default_gpu", "beam_gpu", "antenna"]},

    'eventify' : { 'compiler' : 'clang++',
                 'cmake' : '-DCMAKE_VERBOSE_MAKEFILE=ON -DBACKEND="eventify"',
                 'env' : 'KMP_AFFINITY= ',
                 'prefix' : '',
                 'args' : [ '', '', '', '', ''],
                 'exe_name' : 'minipic.eventify',
                 'threads' : [8, 8, 8, 1, 1],
                 'benchmarks' : [ "default", "beam", "antenna", "Ecst", "Bcst" ]},

    'sycl' : { 'compiler' : 'icpx',
                 'cmake' : '-DCMAKE_VERBOSE_MAKEFILE=ON -DBACKEND="sycl"',
                 'env' : '',
                 'prefix' : '',
                 'args' : [ '', '', '', '', ''],
                 'exe_name' : 'minipic.icpx',
                 'threads' : [8, 8, 8],
                 'benchmarks' : [ "default_gpu", "beam_gpu", "antenna"]},

    'sycl_gpu' : { 'compiler' : 'icpx',
                 'cmake' : '-DCMAKE_VERBOSE_MAKEFILE=ON -DBACKEND="sycl" -DDEVICE="nvidia_v100"',
                 'env' : '',
                 'prefix' : '',
                 'args' : [ '-gpu', '-gpu', '-gpu', '', ''],
                 'exe_name' : 'minipic.icpx',
                 'threads' : [1, 1, 1],
                 'benchmarks' : [ "default_gpu", "beam_gpu", "antenna"]},

    'sycl_pvc' : { 'compiler' : 'icpx',
                 'cmake' : '-DCMAKE_VERBOSE_MAKEFILE=ON -DBACKEND="sycl" -DDEVICE="intel_pvc"',
                 'env' : '',
                 'prefix' : '',
                 'args' : [ '-gpu', '-gpu', '-gpu', '', ''],
                 'exe_name' : 'minipic.icpx',
                 'threads' : [1, 1, 1],
                 'benchmarks' : [ "default_gpu", "beam_gpu", "antenna"]},

    'acpp' : { 'compiler' : 'acpp',
                 'cmake' : '-DCMAKE_VERBOSE_MAKEFILE=ON -DBACKEND="acpp" -DDEVICE="cpu_x86"',
                 'prefix' : '',
                 'args' : [ '', '', '', '', ''],
                 'exe_name' : 'minipic.acpp',
                 'threads' : [8, 8, 8],
                 'benchmarks' : [ "default_gpu", "beam_gpu", "antenna"]},

    'stdpar_gpu' : { 'compiler' : 'nvc++',
                # 'cmake' : '-DCMAKE_VERBOSE_MAKEFILE=ON -DBACKEND="stdpar" -DDEVICE="nvidia_v100"',
                'cmake' : '-DCMAKE_VERBOSE_MAKEFILE=ON -DBACKEND="stdpar" -DDEVICE="nvidia_a100"',
                'env' : '',
                'prefix' : '',
                'args' : [ '', '', '', '', ''],
                'exe_name' : 'minipic.stdpar',
                'threads' : [1, 1, 1],
                'benchmarks' : [ "default_gpu", "beam_gpu", "antenna"]},

    'stdpar_cpu' : { 'compiler' : None,
            'cmake' : '-DCMAKE_VERBOSE_MAKEFILE=ON -DBACKEND="stdpar_cpu"',
            'env' : '',
            'prefix' : '',
            'args' : [ '', '', '', '', ''],
            'exe_name' : 'minipic.stdparcpu',
            'threads' : [8, 8, 8, 1, 1],
            'benchmarks' : [ "default_gpu", "beam", "antenna", "Ecst", "Bcst" ]},
}

config_description = "List of all configurations: \n\n"

for config in configuration_list.keys():
    config_dict = configuration_list[config]
    config_description += "> {}: \n".format(config)
    config_description += "    - compiler:  \n".format(config_dict['compiler'])
    config_description += "    - cmake options:  {}\n".format(config_dict['cmake'])
    config_description += "    - run prefix:  {}\n".format(config_dict['prefix'])

config_description += "    \n\n"

# Command line arguments
parser = argparse.ArgumentParser(prog='Validation process',
                                 description='Custom arguments for the validation process',
                                 epilog=config_description,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-g', '--config', default='openmp', help=' configuration choice: sequential, openmp (default), kokkos, kokkos_gpu, thrust')  
parser.add_argument('-c', '--compiler', help=' custom compiler choice')
parser.add_argument('-b', '--benchmarks', help=' specific benchmark, you can specify several benchmarks with a coma. for instance "default,beam"', default=None)
parser.add_argument('-t', '--threads', help=' default number of threads', default=None)
parser.add_argument('-a', '--arguments', help=' default arguments', default=None)
parser.add_argument('--clean', help=' whether to delete or not the generated files', action='store_true')
parser.add_argument('--no-evaluate', help=' if used, do not evaluate against the reference', action='store_true')
parser.add_argument('--compile-only', help=' if used, only compile the tests', action='store_true')
parser.add_argument('--threshold', help=' threshold for the validation', default=1e-10, type=float)
parser.add_argument('--save-timers', help=' save the timers for each benchmark', action='store_true')
parser.add_argument('--prefix', help=' add custom prefix for the execution, for instance srun', default='')
parser.add_argument('--env', help=' add custom environment variables for the execution, for instance `OMP_PROC_BIND=spread`', default='')
parser.add_argument('--device', help=' select the device type to pass for cmake compilation`', default='')
parser.add_argument('--cmake-args', help=' add custom cmake arguments for the compilation', default='')
parser.add_argument('--backend', help=' select the backend to use', default='')
args = parser.parse_args()

# Selected configuration
configuration = args.config

if (configuration == "omp_task"):
    configuration = "openmp_task"

selected_config = configuration_list[configuration].copy()

# Select compiler
if (args.compiler != None):
    selected_config['compiler'] = args.compiler

clean = args.clean

# Evaluate or not
evaluate = not(args.no_evaluate)

# Compile only
compile_only = args.compile_only

# Save timers
save_timers = args.save_timers

# Select benchmark list
# args.benchmark should have the form of a list of benchmarks separated by a coma like "default,beam"

if (args.benchmarks != None):

    selected_config['benchmarks'] = []

    for b in args.benchmarks.split(","):
        selected_config['benchmarks'].append(b)

# Select threads 
if (args.threads != None):

    # if args.threads is a list of threads

    if ("," in args.threads):

        selected_config['threads'] = []
        for thread in args.threads.split(","):
            selected_config['threads'].append(int(thread))
    
    # if args.threads is a single number
                
    elif (args.threads.isdigit()):

        selected_config['threads'] = []
        for i in range(len(selected_config['benchmarks'])):
            selected_config['threads'].append(int(args.threads))

# Environment
if (args.env != None):
    selected_config['env'] += ' ' + args.env

# Prefix 
if (args.prefix != None):
    selected_config['prefix'] += ' ' + args.prefix

# Select arguments
if (args.arguments != None):

    if ("," in args.arguments):

        selected_config['args'] = []
        for arg in args.arguments.split(","):
            selected_config['args'].append(arg)

    else:

        selected_config['args'] = []
        for i in range(len(selected_config['benchmarks'])):
            selected_config['args'].append(args.arguments)

else:

    selected_config['args'] = []
    for i in range(len(selected_config['benchmarks'])):
        selected_config['args'].append('')

# Select device
if ((args.device != None) and (args.device != '')):

    cmake_args = selected_config['cmake'].split(" ")

    # remove the device option if exists
    cmake_args = [arg for arg in cmake_args if not arg.startswith("-DDEVICE=")]
    cmake_args.append("-DDEVICE={}".format(args.device))
    # add the new device option
    selected_config['cmake'] = " ".join(cmake_args)

# Change backend
if ((args.backend != None) and (args.backend != '')):
    cmake_args = selected_config['cmake'].split(" ")

    # remove the backend option if exists
    cmake_args = [arg for arg in cmake_args if not arg.startswith("-DBACKEND=")]
    cmake_args.append("-DBACKEND={}".format(args.backend))
    # add the new device option
    selected_config['cmake'] = " ".join(cmake_args)

# Add custom cmake arguments
if ((args.cmake_args != None) and (args.cmake_args != '')):
    
    selected_config['cmake'] = selected_config['cmake'] + " " + args.cmake_args

# threshold
threshold = args.threshold

assert(threshold > 0), "Threshold should be positive"

# Check if the configuration is valid:
# - Size of the threads list should be equal or above to the number of benchmarks
assert(len(selected_config['threads']) >= len(selected_config['benchmarks'])), "No enough threads specified ({}) for the number of benchmarks ({})".format(selected_config['threads'], selected_config['benchmarks'])
# - Threads should be positive
assert(all([t > 0 for t in selected_config['threads']])), "Number of threads should be positive ({})".format(selected_config['threads'])


# Get local path
working_dir = os.getcwd() + "/"
root_dir = working_dir + "/../"

# Add validation path to the PYTHONPATH
sys.path.append('{}/validation'.format(root_dir))

# Add root path to the PYTHONPATH
sys.path.append('{}'.format(root_dir))

# Terminal size
try:
    terminal_size = os.get_terminal_size()[0]
except:
    terminal_size = 80

def line(size):
    line = " "
    for i in range(size-1):
         line += "_"
    print(line)

# Get the git hash in variable git_hash
try:
    git_hash = subprocess.check_output(['git', 'rev-parse', 'HEAD'], timeout=60).strip().decode('utf-8')
except:
    git_hash = "unknown"

# Get the git branch in variable git_branch
try:
    git_branch = subprocess.check_output(['git', 'rev-parse', '--abbrev-ref', 'HEAD'], timeout=60).strip().decode('utf-8')
except:
    git_branch = "unknown"

# Get the pipeline id in variable pipeline_id
try:
    pipeline_id = os.environ['CI_PIPELINE_ID']
except:
    pipeline_id = "unknown"

# get branch from detatched head
if git_branch == "HEAD":
    git_branch = "detached"



# Create unique id of the form yyyymmdd_hh_mm_ss_pipeline_id_git_hash
unique_id = time.strftime("%Y%m%d_%H:%M:%S") + "_" + pipeline_id + "_" + git_hash

# timers saving path
timers_path = "/gpfs/workdir/labotm/gitlab-ci_runs/minipic_timers/"

# ________________________________________________________________________________
# Print some info
line(terminal_size)
print(' VALIDATION \n')

print(" Git branch: {}".format(git_branch))
print(" Git hash: {}".format(git_hash))
print(" Pipeline id: {}".format(pipeline_id))
print(" Working dir: {}".format(working_dir))
print(" Selected configuration: {}".format(configuration))
print(" Compiler: {}".format(selected_config['compiler']))
print(" Clean: {}".format(clean))
print(" Save: {}".format(save_timers))
print(" Unique id: {}".format(unique_id))
if compile_only:
    print(" Compile only: {}".format(compile_only))
else:
    print(" Evaluate: {}".format(evaluate))
print(" Threshold: {}".format(threshold))
print(" Prefix: {}".format(selected_config['prefix']))
print(" Env: {}".format(selected_config['env']))
print(" Device: {}".format(args.device))
if args.backend != '' and args.backend != None:
    print(" Backend: {}".format(args.backend))
print(" Cmake args: {}".format(selected_config['cmake']))

# print all benchamrks

print(" Benchmarks: ")
for ib,benchmark in enumerate(selected_config['benchmarks']):
    print("   - {} : {} threads - args: {}".format(benchmark, selected_config['threads'][ib], selected_config['args'][ib] ))

# ________________________________________________________________________________
# Run benchmarks

# Get configuration
compiler = selected_config['compiler']
cmake = selected_config['cmake']
executable_name = selected_config['exe_name']

for ib,benchmark in enumerate(selected_config['benchmarks']):

    env = selected_config['env']
    prefix = selected_config['prefix']
    args = selected_config['args'][ib]

    # bench parameters

    nb_threads = selected_config['threads'][ib]
    bench_dir = working_dir + benchmark + "/"

    print("")
    line(terminal_size)
    print("\n >>> Benchmark `{}` \n".format(benchmark))

    # ____________________________________________________________________________
    # Compilation

    # Create a directory for this benchmark
    if os.path.exists(bench_dir):
        os.system("rm -rf " + bench_dir)
    os.makedirs(bench_dir)

    # Go to the bench directory
    os.chdir(bench_dir)

    # Copy the main from src
    os.system("cp " + root_dir + "src/main.cpp .")

    # Change include
    with open("main.cpp", "r") as file:
	    main_file = file.readlines()

    start_index = 0
    while("Load a setup" not in main_file[start_index]):
        start_index += 1

    # start_index = main_file.index("Load a setup")

    for iline in range(start_index+1, len(main_file), 1):
        if '#include' in main_file[iline]:
            main_file[iline] = '\n'

    main_file[start_index+1] = '#include "{}.hpp" \n'.format(benchmark)

    with open("main.cpp", "w") as file:
        file.writelines(main_file)

    # Copy Cmake file
    # os.system("cp ./../CMakeLists.txt .")

    # Compile

    cmake_command = "cmake ../../ -DCMAKE_BUILD_TYPE=Release -DTEST=ON"
    #cmake_command += " -DCMAKE_INSTALL_PREFIX={}".format(bench_dir)
    if (compiler != None):
        cmake_command += " -DCMAKE_CXX_COMPILER={}".format(compiler)
    cmake_command += " {}".format(cmake)

    print("")
    print("   -> Compilation")
    print("")
    print(cmake_command)

    os.system(cmake_command)
    os.system("make -j 4")

    # Check
    exe_exists = os.path.exists("./{}".format(executable_name))
    assert(exe_exists), "Executable not generated"  

    # ____________________________________________________________________________
    # Execution

    if not(compile_only):

        #os.system("{} ./{} {}".format(prefix, executable_name, args))

        # if benchmark has key threads
        env = " OMP_NUM_THREADS={} ".format(nb_threads) + env

        print("")
        print("   -> Execution ")
        print("")

        subprocess_command = ["{} {} ./{} {}".format(env, prefix, executable_name, args)] # srun numactl --interleave=all 

        print(subprocess_command)

        subprocess.run(subprocess_command, shell=True, check=False)

        # ____________________________________________________________________________
        # Check results

        # Load the corresponding module
        if (os.path.exists(root_dir + "/validation/{}.py".format(benchmark))):

            print("")
            print("   -> Launch the validation process ")
            print("")

            module = importlib.import_module(benchmark, package=None)

            module.validate(evaluate=evaluate, threshold=threshold)

            print("")
            print(" \033[32mBenchmark `{}` tested with success \033[39m".format(benchmark))
            line(terminal_size)

        # ____________________________________________________________________________
        # Read timers (json format)

        if save_timers:

            print("")
            print("   -> Timers")
            print("")

            with open("timers.json") as f:
                raw_timers_dict = json.load(f)

            pic_it_time = raw_timers_dict['final']['main loop (no diags)'][0]

            # Append the following information to the file ci_benchmark_timers.json
            # {
            #   [unique_id] : {
            #      "date" : "yyyy_mm_dd_hh_mm_ss",
            #      "branch" : branch,
            #      "hash" : hash,
            #      "pipeline_id" : pipeline_id,
            #      "times" : {
            #          [iteration] : pic_it_time
            #      }
            # }
            # if the file does not exist, create it

            ci_file_name = timers_path + "/ci_{}_{}_timers.json".format(benchmark, configuration)

            if not os.path.exists(ci_file_name):
                with open(ci_file_name, "w") as f:
                    f.write("{}")
                f.close()

            with open(ci_file_name, "r") as f:
                ci_timers_dict = json.load(f)

            # Read all ci_timers_dict keys and compute an average time:
                
            average_pic_it_time = 0.0
            number_of_entry = len(ci_timers_dict.keys())
            for key in ci_timers_dict.keys():
                average_pic_it_time += ci_timers_dict[key]["times"]["iteration"]

            if number_of_entry > 0:
                average_pic_it_time /= number_of_entry

            print("    Average iteration time: {}".format(average_pic_it_time))
            print("    Current iteration time: {}".format(pic_it_time))

            # Create the new entry

            ci_timers_dict[unique_id] = { "date" : time.strftime("%Y_%m_%d_%H_%M_%S"),
                                        "branch" : git_branch,
                                        "hash" : git_hash,
                                        "pipeline_id" : pipeline_id,
                                        "times" : { "iteration" : pic_it_time } }
            
            # Save the new entry

            with open(ci_file_name, "w") as f:
                json.dump(ci_timers_dict, f, indent=4)
            
            f.close()

            print("    Timers saved in {}".format(ci_file_name))
            
        
                       
        


    # ____________________________________________________________________________

    if clean and os.path.exists(bench_dir):

        print("")
        print("   -> Cleaning")

        os.system("rm -rf " + bench_dir)
