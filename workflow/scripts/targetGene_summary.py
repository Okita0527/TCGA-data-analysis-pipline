import pandas as pd
import click
#python targetGene_summary.py -r results/TCGA-LUAD.raw.count.tn.xls -n results/TCGA-LUAD.normalization.count.tn.xls -c results/TCGA-LUAD_clinic.xls -t ENSG00000146648 -o results/EGFR_summary.xls


def get_targetgene_index(count_data,target_gene):
    for i in range(0, len(count_data["Unnamed: 0"])):
        id_v = count_data.loc[i]["Unnamed: 0"]
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

def get_count(clinic_data,count_data):
    tumor_count = {}
    normal_count = {}
    drop_list = []
    for i in range(0,len(clinic_data.index)):
        patient_barcode = clinic_data.loc[i]["submitter_id"]
        flags = 0
        for all_barcode in count_data.index[1:]:
            barcode_3 = get_barcode3(all_barcode)
            if "-".join(barcode_3.split("-")[:-1]) == patient_barcode:
                if "tumor" in barcode_3:
                    tumor_count[patient_barcode] = count_data[all_barcode]
                    flags +=1
                elif "normal" in barcode_3:
                    normal_count[patient_barcode] = count_data[all_barcode]
                    flags += 1
                else:
                    print(f"{all_barcode} is not belong to tumor tissue or normal tissue, please check it")
            if flags == 2:
                break

        if flags < 2:
            print(f"The count corresponding to {patient_barcode} does not exist ,index '{i}' will be deleted")
            drop_list.append(i)


    #output_dataframe = clinic_data.drop(drop_list, axis=0)
    #print(len(output_dataframe.index))
    return tumor_count,normal_count,drop_list

def merge_dataframe(clinic_data,drop_list,raw_tumor_count,raw_normal_count,normalization_tumor_count,normalization_normal_count):
    output_DataFrame = clinic_data.drop(drop_list, axis=0)
    output_DataFrame = output_DataFrame.reset_index(drop=True)
    output_DataFrame.insert(loc=2, column="Ratio of Normalized count(Cancer/Normal)", value="NA")
    output_DataFrame.insert(loc=3, column="Raw Count in Cancer", value="NA")
    output_DataFrame.insert(loc=4, column="Raw Count in Normal", value="NA")
    output_DataFrame.insert(loc=5, column="Normalized Count in Cancer", value="NA")
    output_DataFrame.insert(loc=6, column="Normalized Count in Normal", value="NA")
    output_DataFrame.insert(loc=7, column="Rank", value="NA")

    for i in output_DataFrame.index:
        patient_barcode = output_DataFrame.loc[i,"submitter_id"]
        output_DataFrame.loc[i,"Raw Count in Cancer"] = raw_tumor_count[patient_barcode]
        output_DataFrame.loc[i,"Raw Count in Normal"] = raw_normal_count[patient_barcode]
        output_DataFrame.loc[i,"Normalized Count in Cancer"] = normalization_tumor_count[patient_barcode]
        output_DataFrame.loc[i,"Normalized Count in Normal"] = normalization_normal_count[patient_barcode]
        output_DataFrame.loc[i,"Ratio of Normalized count(Cancer/Normal)"] = \
        float(normalization_tumor_count[patient_barcode]) / float(normalization_normal_count[patient_barcode])
    output_DataFrame['Rank'] = output_DataFrame["Ratio of Normalized count(Cancer/Normal)"].rank()
    return output_DataFrame

#python targetGene_summary.py -r results/TCGA-LUAD.raw.count.tn.xls -n results/TCGA-LUAD.normalization.count.tn.xls -c results/TCGA-LUAD_clinic.xls -t NSG00000146648 -o results/TargetGene_summary.xls

@click.command()
@click.option('-r', 'raw_count_file', type=str, default="./results/TCGA-LUAD.raw.count.tn.xls", show_default=True,
              help='输入的count文件，由“extract_tn.py”提供')
@click.option('-n', 'normalization_count_file', type=str, default="./results/TCGA-LUAD.normalization.count.tn.xls", show_default=True,
              help='输入的count文件，由“extract_tn.py”提供')
@click.option('-c', 'clinic_file', type=str, default="./results/TCGA-LUAD_clinic.xls", show_default=True,
              help='输入的count文件，由“download_count_clinic.R”提供')
@click.option('-t', 'target_gene', type=str, default="ENSG00000146648", help='目的基因，使用ensembl的ID')
@click.option('-o', 'output_file', type=str, default="./results/EGFR_summary.xls", show_default=True,
              help='靶基因汇总文件，包括靶基因所对应样本的count数及临床信息')

def main(raw_count_file,normalization_count_file,clinic_file,target_gene,output_file):
    raw_count = pd.read_csv(raw_count_file, sep="\t")
    normalization_count = pd.read_csv(normalization_count_file, sep="\t")
    clinic_data = pd.read_csv(clinic_file, sep="\t")

    targetGene_index = get_targetgene_index(raw_count,target_gene)
    targetGene_raw_count = raw_count.loc[targetGene_index]  # index is barcode,label is target_gene
    targetGene_normalization_count = normalization_count.loc[targetGene_index]  # index is barcode,label is target_gene
    raw_tumor_count,raw_normal_count,drop_list = get_count(clinic_data,targetGene_raw_count)  # dictionary,such as {'TCGA-55-6986': 4077}
    normalization_tumor_count,normalization_normal_count,drop_list = get_count(clinic_data,targetGene_normalization_count)

    output_DataFrame = merge_dataframe(clinic_data, drop_list,
                                       raw_tumor_count, raw_normal_count,
                                       normalization_tumor_count, normalization_normal_count)



    output_DataFrame = output_DataFrame.drop(["project","bcr_patient_barcode"], axis=1)
    output_DataFrame.to_csv(output_file, sep="\t",index=False)


if __name__ == "__main__":
    main()





