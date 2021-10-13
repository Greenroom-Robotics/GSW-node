# A script to read the methods from the c++ header files. 

import readline

def input_with_prefill(prompt, text):
    def hook():
        readline.insert_text(text)
        readline.redisplay()
    readline.set_pre_input_hook(hook)
    result = input(prompt)
    readline.set_pre_input_hook()
    return result

file_input = input_with_prefill("File path:", "")

# Get a list of methods from a header file for use in the bindings
_keyword = "gsw_"

print(file_input)
with open(file_input) as fp:
    for line in fp:
        if _keyword in line:
            start_index = line.index(_keyword)
            remove_before = line[start_index:len(line)]
            method = remove_before.split("(")[0]
            print(method)
