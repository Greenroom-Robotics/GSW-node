# This script is used to fix an import error in the generated bindings. 
# Do not run this script manually, it is used during the install process. 
# This is the outcome of using a forked version of nbind to allow use with
# a newer version of node. 

print("Fixing nbind import in lib-types.d.ts")

file_name = "./lib-types.d.ts"
file_content = ""

with open(file_name, "r") as f:
    first_row = f.readline()
    replaced_import = first_row.replace("nbind", "@mcesystems/nbind")

    file_content += replaced_import
    for line in f:
        file_content += line

with open(file_name, "w") as f:
    f.write(file_content)


