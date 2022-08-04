import pandas as pd
import click
###usage: python extract_tn.py -i results/TCGA-LUAD.raw.count.csv -t ENSG00000146648 -o results/TCGA-LUAD.raw.count.tn.xls



def get_targetgene_index(count_data,target_gene):
    for i in range(0, len(count_data["Unnamed: 0"])):
        id_v = count_data.iloc[i]["Unnamed: 0"]
        id = id_v.split(".")[0]
        if id == target_gene:
            return (i)

def get_barcode3(all_barcode):
    tumor_types = ["01", "02", "03", "04", "05", "06", "07", "08", "09"]
    normal_types = ["10", "11", "12", "13", "14", "15", "16", "17", "18", "19"]
    control_samples = ["20", "21", "22", "23", "24", "25", "26", "27", "28", "29"]
    barcode3 = all_barcode.split("-")[:3]
    barcode3 = "-".join(barcode3)
    if all_barcode.split("-")[3][:2] in tumor_types:
        barcode3 += "-tumor"
    elif all_barcode.split("-")[3][:2] in normal_types:
        barcode3 += "-normal"
    elif all_barcode.split("-")[3][:2] in control_samples:
        barcode3 += "-control"
    else:
        print(f"Sample label in {all_barcode} does not conform to naming conventions which made by GDC, please check it")
    return barcode3


def remove_duplicate_barcode(count_data,targetGene_index):
    new_barcode = {}
    for each_barcode in count_data.columns[1:]:
        barcode3 = get_barcode3(each_barcode)
        if not barcode3 in new_barcode.keys():
            new_barcode[barcode3] = each_barcode
        else:
            old_count = count_data.iloc[targetGene_index][new_barcode[barcode3]]
            new_count = count_data.iloc[targetGene_index][each_barcode]
            if int(new_count) > int(old_count):
                new_barcode[barcode3] = each_barcode
                print(f"------------------------------------------------------------\n"
                      f"target gene count of {new_barcode[barcode3]} is {old_count}\n"
                      f"target gene count of {each_barcode} is {new_count}\n"
                      f"{new_barcode[barcode3]} has been delete")
            else:
                print(f"------------------------------------------------------------\n"
                      f"target gene count of {new_barcode[barcode3]} is {old_count}\n;"
                      f"target gene count of {each_barcode} is {new_count}\n"
                      f"{each_barcode} has been delete")
    return new_barcode

def find_Samplepair(barcode_dict):
    new_barcode = []
    barcode_3 = []
    for i in barcode_dict:
        i = i.split("-")[:-1]
        ii = "-".join(i)
        barcode_3.append(ii)

    for j in barcode_dict:
        jj = j.split("-")[:-1]
        jj ="-".join(jj)
        numbers = barcode_3.count(jj)
        if numbers >= 2:
            new_barcode.append(barcode_dict[j])
    return new_barcode


@click.command()
@click.option('-i', 'input_file', type=str, default="results/TCGA-LUAD.raw.count.csv", show_default=True,
              help='输入的count文件，由上一步“download_count_clinic.R”提供')
@click.option('-t', 'target_gene', type=str, default="ENSG00000146648", help='目的基因，使用ensembl的ID')
@click.option('-o', 'output_file', type=str, default="results/TCGA-LUAD.raw.count.tm.xls", show_default=True,
              help='输出的count文件，仅保留成对样本（癌与癌旁）且在ID重复时删除靶基因对应count较小的样本')

def main(input_file,target_gene,output_file):
    count_data = pd.read_csv(input_file, sep=",")
    targetGene_index = get_targetgene_index(count_data,target_gene)
    barcode_del_duplicate = remove_duplicate_barcode(count_data,targetGene_index)
    new_barcode = find_Samplepair(barcode_del_duplicate)
    tn_count_data = count_data[new_barcode]
    tn_count_data.index = count_data["Unnamed: 0"]
    tn_count_data.index.name = ''
    tn_count_data.to_csv(output_file, sep="\t")

if __name__ == "__main__":
    main()