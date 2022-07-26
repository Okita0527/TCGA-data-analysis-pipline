def get_barcode(input_filename):
    list1 = []
    input_file = open(input_filename,"r").readlines()
    for eachline in input_file:
        eachline = eachline.strip()
        list1.append(eachline)
    return list1

old_barcode = get_barcode("old_barcode.txt")

def dataset_to_xls(dataset,output_filename):
    output_file = open(output_filename,"w")
    for eachline in dataset:
        eachline = list(map(lambda x: str(x),eachline))
        eachline = "\t".join(eachline)
        output_file.write("{}\n".format(eachline))
    output_file.close()
def read_xls(input_filename):
    input_file = open(input_filename, "r").readlines()
    raw_dataset = list(map(lambda x : x.strip().split("\t"),input_file))
    return raw_dataset

raw_dataset = read_xls("fpkm_unstranded.xls")
raw_headline = raw_dataset[0]
raw_headline.insert(1," ")

locations = []
for i in old_barcode:
    for j in range(0,len(raw_headline)):
        jj = raw_headline[j]
        #jj = "-".join(jj.split("-")[:3])
        if i == jj:
            locations.append(j)
            break

new_dataset = []
for eachline in raw_dataset:
    new_line =[]
    for each_location in locations:
        new_line.append(eachline[each_location])
    new_dataset.append(new_line)

dataset_to_xls(new_dataset,"fpkm.xls")
