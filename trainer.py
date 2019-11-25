from gor import GOR3

import csv

gortrain = GOR3()
with open("2018-06-06-ss.cleaned.csv") as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=",")
    line_count = 0
    for row in csv_reader:
        if line_count == 0:
            print(f'Column names are {", ".join(row)}')
            line_count += 1
        else:
            if(row[6].strip() == "False"):
                # print("sequence = " + row[2] + ", structure = " + row[4])
                gortrain.train(row[2], row[4])
            line_count += 1
        # if line_count >= 100:
        #     break
    print(f'Processed {line_count} lines.')
gortrain.save_model()