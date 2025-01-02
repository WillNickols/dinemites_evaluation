from anadama2 import Workflow
import os
import itertools
import copy

workflow = Workflow(version="0.1", description="Evaluate novel infection methodology")
workflow.add_argument("cores", desc="The number of CPU cores allocated to the job", type=int, default=1)
workflow.add_argument("mem", desc="The memory in megabytes allocated to run the command", type=int, default=10000)
workflow.add_argument("time", desc="The time in minutes allocated to run the command", type=int, default=120)
workflow.add_argument('evaluation', desc="Which evaluation to run", default="tmp")
args = workflow.parse_args()
this_directory = str(os.path.dirname(os.path.realpath(__file__))).rstrip('/') + '/'

# grid
memory = args.mem
cores = args.cores
partition = args.grid_partition
time = args.time
evaluation = args.evaluation

output = os.path.abspath(args.output.rstrip("/")) + '/'
os.makedirs(output, exist_ok=True)
nIterations = 50


# import parameters from
if evaluation == 'qpcr':
    parameters_file = 'scripts/config/qpcr_testing.txt'
elif evaluation == 'locus':
    parameters_file = 'scripts/config/locus_testing.txt'
elif evaluation == 'other':
    parameters_file = 'scripts/config/other_testing.txt'
else:
    raise ValueError("evaluation not valid")

run_only_args = ['model', 'multiple_imputation']

with open(parameters_file, 'r') as file:
    lines = file.readlines()

def parse_parameters(lines):
    param_dict = {}
    
    for line in lines:
        # Split by ':'
        key, value = line.split(':', 1)
        key = key.strip()
        value = value.strip().split(':')
        
        # Get the combining method and parameters
        combine_method = value[0].strip()  # 'one' or 'all'
        options = value[1].strip().split()  # Options after the last colon
        
        # Store the options in the dictionary
        param_dict[key] = {
            'method': combine_method,
            'options': options
        }
    
    return param_dict

def generate_combinations(param_dict):
    combinations = []
    
    # Prepare lists for options
    one_options = {}
    all_options = {}

    for key, params in param_dict.items():
        method = params['method']
        options = params['options']
        
        if method == 'one':
            # Use the first option by default
            one_options[key] = options[0]  # Only take the first option
        elif method == 'all':
            # Use all options
            all_options[key] = options
    
    # If there are no 'all' parameters, we just return one combination
    if not all_options:
        combinations.append(one_options)
        return combinations

    # Generate combinations by iterating through all 'all' parameters
    all_keys = list(all_options.keys())
    all_values = [all_options[k] for k in all_keys]

    # For each combination of all parameters
    for prod in itertools.product(*all_values):
        for key, first_value in one_options.items():
            # Hold other "one" parameters constant and vary the current one
            current_combination = one_options.copy()  # Start with 'one' options
            current_combination[key] = first_value  # Hold the first option constant
            
            for option in param_dict[key]['options']:
                current_combination[key] = option  # Change to the current option

                # Add combinations of 'all' parameters
                for i, all_key in enumerate(all_keys):
                    current_combination[all_key] = prod[i]  # Add the combinations of 'all' parameters

                combinations.append(current_combination.copy())

    return combinations

def remove_repeated_values(dict_list):
    seen = set()
    unique_dicts = []
    duplicated_dicts = []

    for d in dict_list:
        values = tuple(d.items())  # Convert dictionary values to a tuple
        if values not in seen:  # Check if the values have been seen before
            seen.add(values)  # Add the tuple to the set if it's new
            unique_dicts.append(d)  # Add the dictionary to the result list
        else:
            duplicated_dicts.append(d)

    if evaluation == 'high':
        for d in duplicated_dicts:
            for u in unique_dicts:
                if tuple(u.items()) == tuple(d.items()):
                    unique_dicts.remove(u)


    return unique_dicts

# Process input data
param_dict = parse_parameters(lines)
combinations = generate_combinations(param_dict)

combinations = remove_repeated_values(combinations)

ordered_keys = [line.split(':', 1)[0] for line in lines]

input_files_generated = set()
for param in combinations:
    for iteration in range(1, nIterations + 1):
        input_file = output + 'input/' + '_'.join(str(param[key]) for key in ordered_keys if key not in run_only_args) + f'_{iteration}.tsv'
        
        new_command = 'Rscript scripts/all_model_data_generation.R' + \
            ' ' + ' '.join([f'--{key} {param[key]}' for key in ordered_keys if key not in run_only_args]) + \
            f' --iteration {iteration}' + \
            f' --output {output}'

        if not os.path.exists(input_file) and input_file not in input_files_generated:
            workflow.add_task_gridable(actions=new_command,
                targets=[input_file],
                time=10,
                mem=2000,
                cores=1,
                partition=partition
            )
        
        input_files_generated.add(input_file)

        output_file = output + 'output/' + '_'.join(str(param[key]) for key in ordered_keys) + f'_{iteration}.tsv'

        new_command = 'Rscript scripts/all_model_testing.R' + \
            ' ' + ' '.join([f'--{key} {param[key]}' for key in ordered_keys]) + \
            f' --iteration {iteration}' + \
            f' --output {output}'

        if not os.path.exists(output_file):
            workflow.add_task_gridable(actions=new_command,
            depends = [input_file],
            targets = [output_file],
            time=time,
            mem=memory,
            cores=cores,
            partition=partition
            )

workflow.go()


