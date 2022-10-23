import csv
import os


files = [data for data in os.walk(os.path.dirname(__file__) + "/smiles")][0][2]

for f in files:
    base_name = f.split(".")[0]
    text_file_name = os.path.dirname(__file__) + "/smiles/" + str(f)
    csv_file_name = os.path.dirname(__file__) + "/smiles/" + base_name + ".csv"

    with open(csv_file_name, "w+", newline="") as csv_file_target:
        csv_file = csv.writer(csv_file_target, delimiter=",")
        csv_file.writerow(["smiles", "refcode"])

        with open(text_file_name, "r+") as text_file_to_csv_convert:

            for line in text_file_to_csv_convert.readlines():
                line_data = line.strip().split()
                if base_name == "smiles":
                    csv_file.writerow(line_data[1:])
                else:
                    csv_file.writerow(line_data)
