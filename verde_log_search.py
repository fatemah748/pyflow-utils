import sys
import shutil
import os

def main(argv = None):
    if argv is None:
       argv = sys.argv

    if len(argv) != 3:
       sys.exit("Usage: verde_log_search.py path/to/search/directory /path/to/copy/directory")

    search_dir = argv[1]
    copy_dir = argv[2]

    exclude = ["ts", "pm7", "rot", "sp-dft", "input", "conf"]
    if os.path.isdir(copy_dir):
        print(copy_dir + " exists")
    else:
        os.mkdir(copy_dir)


    print("copying all log files from " + search_dir + " to " + copy_dir)
    print("excluding files that contain: " + str(exclude))

    # walking through all directories and subdirectories
    for root, dirs, files in os.walk(search_dir):
        for file in files:
            if file.endswith('.log'):
                # creates full file path
                copy_file = root + "/" + file
                # searches file string for excluded strings
                test = [s for s in exclude if (s in file)]
                if test == []:
                    # checks if the file terminated normally by reading last line
                    if 'Normal termination' in open(copy_file).readlines()[-1]:
                        print("copied file: " + file)
                        # copies file to directory
                        shutil.copy(copy_file, copy_dir)
                    else:
                        print("file did not terminate normally: " + file)
                else:
                    print("excluded file: " + file)

    return 0

if __name__ == '__main__':
   main()