import pandas as pd
import sklearn
from sklearn import metrics
import matplotlib.pyplot as plt

import glob
import os
import shutil
import time
import subprocess

FOLDER = "SIMULATED_SEQ"

def run_guidance_on_set(folder):
    times = []
    count = 0
    for fasta_path in glob.glob(folder + "/*.fas"):
        fasta_file = fasta_path.split("/")[1]
        code = fasta_file.split(".fas")[0]
        try:
            shutil.copy(os.path.join(os.getcwd(), folder, fasta_file), os.path.join(os.getcwd(), fasta_file))
            script = os.path.join(os.getcwd(), "guidance_main.py")
            # script_perl = os.path.join("/Users/kpolonsky/Documents/GUIDANCE-guidance.v2.02/www/Guidance/", "guidance.pl")
            # cmd = f"python3 {script} --seqFile  /Users/kpolonsky/PycharmProjects/HoT_Py/{fasta_path} --msaProgram MAFFT --seqType aa --outDir /Users/kpolonsky/PycharmProjects/HoT_Py/{code}/ --proc_num 8"
            cmd = f"python3 {script} --seqFile  /Users/kpolonsky/PycharmProjects/HoT_Py/{fasta_path} --msaProgram MAFFT --seqType aa --outDir /Users/kpolonsky/PycharmProjects/HoT_Py/{code}/"
            # cmd = f"perl {script_perl} --seqFile  /Users/kpolonsky/PycharmProjects/HoT_Py/{fasta_path} --msaProgram MAFFT --seqType aa --outDir /Users/kpolonsky/PycharmProjects/HoT_Py/{code}/"
            # Start timer
            start_time = time.time()
            subprocess.run(cmd, shell=True, check=True)
            # End timer
            end_time = time.time()
            time_taken = end_time - start_time
            times.append(time_taken)
            count += 1
        except Exception as e:
            print(f"failed to run {fasta_file} : {e}\n")
        os.system(f"rm {os.path.join(os.getcwd(), fasta_file)}")
    print(f"Average time per run of guidance analysis: \n {sum(times)/count/60:.2f} minutes")
    print(f"Times taken for guidance analysis: \n {times}")
    n, bins, patches = plt.hist(times)
    plt.show()


def run_msa_set_score_on_set():
    # for fasta_path in glob
    # for fasta_path in glob.glob("/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/SIMULATED_SEQ" + "/*.fas"):
    for fasta_path in glob.glob("/Users/kpolonsky/PycharmProjects/HoT_Py/SIMULATED_SEQ" + "/*.fas"):
        fasta_file = fasta_path.split("/")[6]
        code = fasta_file.split(".fas")[0]
        # MSA_to_score = f"/Users/kpolonsky/PycharmProjects/HoT_Py/{code}/MSA.MAFFT.aln.With_Names"
        MSA_to_score = f"/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/ORIGINAL_PERL_SEMPHY_115_1_CPU/{code}/MSA.MAFFT.aln.With_Names"
        true_MSA = f"/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/SIMULATED_MSA/{code}_TRUE.fas"
        output_files_prefix = "scored_vs_true"
        try:
            # ./msa_set_score  file_with_MSA_to_score   output_files_prefix   -m file_with_a_single_alternative_MSA
            cmd = f"./programs/msa_set_score/msa_set_score   {MSA_to_score}   {output_files_prefix}   -m {true_MSA}"
            print(cmd + "\n")
            pr = open("/Users/kpolonsky/PycharmProjects/HoT_Py/" + code + "_setScoreLog.txt", "w")
            p = subprocess.Popen(cmd, stdout=pr, shell=True)
            p.wait()
            os.makedirs(os.path.dirname(
                f"/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/TRUE_SCORE_RESULTS_KSENIA_115msa_PERL_SEMPHY_1_CPU/Scored_{code}/"),
                        exist_ok=True)
            for file in glob.glob('/Users/kpolonsky/PycharmProjects/HoT_Py/scored_vs_true*'):
                shutil.move(file,
                            f"/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/TRUE_SCORE_RESULTS_KSENIA_115msa_PERL_SEMPHY_1_CPU/Scored_{code}/")
            # subprocess.run(cmd, shell=True, check=True)
        except Exception as e:
            print(f"failed to run {fasta_file} : {e}\n")


def calculateAuc(predictedfile, truefile, code):
    df1 = pd.read_csv(truefile, sep='\s+', header=0, names=['coln', 'rn1', 'rn2', 'truescore'], low_memory=False)
    df2 = pd.read_csv(predictedfile, sep='\s+', header=0, names=['coln', 'rn1', 'rn2', 'predcore'], low_memory=False)
    result = pd.merge(df1, df2, on=['coln', 'rn1', 'rn2'], how="left")
    result['predcore'] = result['predcore'].fillna(0)
    result = result.dropna()
    y_true = result["truescore"].to_numpy()
    y_score = result["predcore"].to_numpy()
    fpr, tpr, thresholds = sklearn.metrics.roc_curve(y_true, y_score, pos_label=None, sample_weight=None,
                                                     drop_intermediate=True)
    roc_auc = metrics.auc(fpr, tpr)

    # plt.title('Receiver Operating Characteristic')
    # plt.plot(fpr, tpr, 'b', label='AUC = %0.2f' % roc_auc)
    # plt.legend(loc='lower right')
    # plt.plot([0, 1], [0, 1], 'r--')
    # plt.xlim([0, 1])
    # plt.ylim([0, 1])
    # plt.ylabel('True Positive Rate')
    # plt.xlabel('False Positive Rate')
    # # plt.show()
    # plt.savefig(f'/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/ROCs/{code}.png')
    return roc_auc


def calculate_average_auc():
    auc = []
    for fasta_path in glob.glob("/Users/kpolonsky/PycharmProjects/HoT_Py/SIMULATED_SEQ" + "/*.fas"):
    # for fasta_path in glob.glob("/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/SIMULATED_SEQ" + "/*.fas"):
        fasta_file = fasta_path.split("/")[6]
        code = fasta_file.split(".fas")[0]
        predicted_file = f"/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/SEMPHY_PYTHON_113_8_CPUs/{code}/MSA.MAFFT.Guidance2_res_pair.scr"
        true_file = f"/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/TRUE_SCORE_RESULTS_PYTHON_113msa_SEMPHY_8_CPUs/Scored_{code}/scored_vs_true_res_pair.scr"
        roc_auc = calculateAuc(predicted_file, true_file, code)
        auc.append(roc_auc)

    # print(" ".join(str(auc)) + "\n")
    print(auc)
    print(f" Average AUC is: {sum(auc) / len(auc)}\n")
    n, bins, patches = plt.hist(auc)
    plt.show()


def calculateAuc_columns(predictedfile, truefile, code):
    df1 = pd.read_csv(truefile, sep='\s+', header=0, names=['col_number', 'col_score_true'], low_memory=False)
    df2 = pd.read_csv(predictedfile, sep='\s+', header=0, names=['col_number', 'col_score_pred'], low_memory=False)
    result = pd.merge(df1, df2, on=['col_number'], how="left")
    result['col_score_pred'] = result['col_score_pred'].fillna(0)
    result = result.dropna()
    y_true = result["col_score_true"].to_numpy()
    y_score = result["col_score_pred"].to_numpy()
    fpr, tpr, thresholds = sklearn.metrics.roc_curve(y_true, y_score, pos_label=None, sample_weight=None,
                                                     drop_intermediate=True)
    roc_auc = metrics.auc(fpr, tpr)

    # plt.title('Receiver Operating Characteristic')
    # plt.plot(fpr, tpr, 'b', label='AUC = %0.2f' % roc_auc)
    # plt.legend(loc='lower right')
    # plt.plot([0, 1], [0, 1], 'r--')
    # plt.xlim([0, 1])
    # plt.ylim([0, 1])
    # plt.ylabel('True Positive Rate')
    # plt.xlabel('False Positive Rate')
    # # plt.show()
    # plt.savefig(f'/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/ROCs/{code}.png')
    return roc_auc


def calculate_average_auc_column_score():
    auc = []
    for fasta_path in glob.glob("/Users/kpolonsky/PycharmProjects/HoT_Py/SIMULATED_SEQ" + "/*.fas"):
    # for fasta_path in glob.glob("/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/SIMULATED_SEQ" + "/*.fas"):
        fasta_file = fasta_path.split("/")[6]
        code = fasta_file.split(".fas")[0]
        predicted_file = f"/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/IQTREE_PYTHON_113_8_CPUs_N0/{code}/MSA.MAFFT.Guidance2_col_col.scr"
        true_file = f"/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/TRUE_SCORE_RESULTS_PYTHON_113msa_IQTREE_8_CPUs_N0/Scored_{code}/scored_vs_true_col_col.scr"
        roc_auc = calculateAuc_columns(predicted_file, true_file, code)
        auc.append(roc_auc)

    # print(" ".join(str(auc)) + "\n")
    print(auc)
    print(f" Average AUC is: {sum(auc) / len(auc)}\n")
    n, bins, patches = plt.hist(auc)
    plt.show()


def save_line_as_tree(line, destination_folder, code):
    with open(os.path.join(destination_folder, f"{code}_indelible_true_tree.ctrl"), "w") as file:
        parts = line.split('TREE')
        tree = parts[2].lstrip()
        file.write(tree)

def get_indelible_true_tree(source_path, destination_folder, code):
    with open(source_path, "r") as source_file:
        # lines = source_file.readlines()
        # for line in lines:
        #     if line.startswith("[TREE]"):
        while True:
            line = source_file.readline()
            if "[TREE]" not in line:
                continue
            else:
                save_line_as_tree(line, destination_folder, code)
                return

def extract_rf(file):
    with open(file, "r") as source_file:
        header = source_file.readline()
        line = source_file.readline()
        data = line.split()
        rf = data[1]
        return rf

def prepare_directories_for_rf():
    # df = pd.DataFrame(columns=["rf1", "rf2", "rf3"])
    for fasta_path in glob.glob("/Users/kpolonsky/PycharmProjects/HoT_Py/SIMULATED_SEQ" + "/*.fas"):
        fasta_file = fasta_path.split("/")[6]
        code = fasta_file.split(".fas")[0]

        # createdirectories
        os.makedirs(os.path.dirname(
            f"/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/RF_DISTANCE/{code}/"),
            exist_ok=True)
        os.makedirs(os.path.dirname(
            f"/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/RF_DISTANCE/{code}/MAFFT/"),
            exist_ok=True)
        os.makedirs(os.path.dirname(
            f"/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/RF_DISTANCE/{code}/Indellible_simulated_alignment/"),
            exist_ok=True)
        os.makedirs(os.path.dirname(
            f"/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/RF_DISTANCE/{code}/MAFFT_after_column_deletions/"),
            exist_ok=True)
        os.makedirs(os.path.dirname(
            f"/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/RF_DISTANCE/{code}/True_tree/"),
            exist_ok=True)
        os.makedirs(os.path.dirname(
            f"/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/RF_DISTANCE/{code}/RFs/"),
            exist_ok=True)

        # copy MSAs
        msa_file_path = f"/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/IQTREE_PYTHON_113_8_CPUs_N0/{code}/MSA.MAFFT.aln.With_Names"
        shutil.copy(msa_file_path, f"/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/RF_DISTANCE/{code}/MAFFT/")
        msa_wo_columns_file_path = f"/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/IQTREE_PYTHON_113_8_CPUs_N0/{code}/MSA.MAFFT.Without_low_SP_Col.With_Names"
        shutil.copy(msa_wo_columns_file_path,
                    f"/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/RF_DISTANCE/{code}/MAFFT_after_column_deletions/")
        indelible_msa_path = f"/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/SIMULATED_MSA/{code}_TRUE.fas"
        shutil.copy(indelible_msa_path,
                    f"/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/RF_DISTANCE/{code}/Indellible_simulated_alignment/")

        # get true tree file out of ctrl
        indelible_tree_source = f"/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/indelible_ctrls/{code}_indelible.ctrl"
        destination_folder = f"/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/RF_DISTANCE/{code}/True_tree/"
        get_indelible_true_tree(indelible_tree_source, destination_folder, code)

        # run IQtree
        cmd = f"cd /Users/kpolonsky/Documents/iqtree/bin; ./iqtree2 -s /Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/RF_DISTANCE/{code}/MAFFT/MSA.MAFFT.aln.With_Names"
        subprocess.run(cmd, shell=True, check=True)
        cmd = f"cd /Users/kpolonsky/Documents/iqtree/bin; ./iqtree2 -s /Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/RF_DISTANCE/{code}/MAFFT_after_column_deletions/MSA.MAFFT.Without_low_SP_Col.With_Names"
        subprocess.run(cmd, shell=True, check=True)
        cmd = f"cd /Users/kpolonsky/Documents/iqtree/bin; ./iqtree2 -s /Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/RF_DISTANCE/{code}/Indellible_simulated_alignment/{code}_TRUE.fas"
        subprocess.run(cmd, shell=True, check=True)

        # run IQtree to calculate RFs
        # rf_file = f"/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/RF_DISTANCE/{code}/True_tree/{code}_indelible_true_tree.ctrl.rfdist"
        # rf_log = f"/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/RF_DISTANCE/{code}/True_tree/{code}_indelible_true_tree.ctrl.log"

    #     # regular MAFFT MSA
    #     cmd = f"cd /Users/kpolonsky/Documents/iqtree/bin; ./iqtree2 -rf /Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/RF_DISTANCE/{code}/MAFFT/MSA.MAFFT.aln.With_Names.treefile /Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/RF_DISTANCE/{code}/True_tree/{code}_indelible_true_tree.ctrl"
    #     subprocess.run(cmd, shell=True, check=True)
    #     shutil.move(rf_file,
    #                 f"/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/RF_DISTANCE/{code}/RFs/{code}_MAFFT_vs_true.rfdist")
    #     shutil.move(rf_log,
    #                 f"/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/RF_DISTANCE/{code}/RFs/{code}_MAFFT_vs_true.rfdist.log")
    #
    #     # MSA without columns
    #     cmd = f"cd /Users/kpolonsky/Documents/iqtree/bin; ./iqtree2 -rf /Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/RF_DISTANCE/{code}/MAFFT_after_column_deletions/MSA.MAFFT.Without_low_SP_Col.With_Names.treefile /Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/RF_DISTANCE/{code}/True_tree/{code}_indelible_true_tree.ctrl"
    #     subprocess.run(cmd, shell=True, check=True)
    #     shutil.move(rf_file,
    #                 f"/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/RF_DISTANCE/{code}/RFs/{code}_MAFFT_without_SP_Col_vs_true.rfdist")
    #     shutil.move(rf_log,
    #                 f"/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/RF_DISTANCE/{code}/RFs/{code}_MAFFT_without_SP_Col_vs_true.rfdist.log")
    #
    #     # simulated MSA
    #     cmd = f"cd /Users/kpolonsky/Documents/iqtree/bin; ./iqtree2 -rf /Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/RF_DISTANCE/{code}/Indellible_simulated_alignment/{code}_TRUE.fas.treefile /Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/RF_DISTANCE/{code}/True_tree/{code}_indelible_true_tree.ctrl"
    #     subprocess.run(cmd, shell=True, check=True)
    #     shutil.move(rf_file,
    #                 f"/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/RF_DISTANCE/{code}/RFs/{code}_indel_simulated_vs_true.rfdist")
    #     shutil.move(rf_log,
    #                 f"/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/RF_DISTANCE/{code}/RFs/{code}_indel_simulated_vs_true.rfdist.log")
    #
    #     # extract RFs from the results
    #     rf1 = extract_rf(f"/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/RF_DISTANCE/{code}/RFs/{code}_MAFFT_vs_true.rfdist")
    #     rf2 = extract_rf(f"/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/RF_DISTANCE/{code}/RFs/{code}_MAFFT_without_SP_Col_vs_true.rfdist")
    #     rf3 = extract_rf(f"/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/RF_DISTANCE/{code}/RFs/{code}_indel_simulated_vs_true.rfdist")
    #
    #     row = {"rf1": [rf1], "rf2": [rf2], "rf3": [rf3]}
    #     df2 = pd.DataFrame(row)
    #     df = pd.concat([df, df2], ignore_index=True)
    # df.to_csv("/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/RF_DISTANCE/rfs.csv", index=False)
    # print(df)

# def test():
#     df = pd.DataFrame(columns=["rf1", "rf2", "rf3"])
#     code= 'CLDN1'
#     # extract RFs from the results
#     rf1 = extract_rf(
#         f"/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/RF_DISTANCE/{code}/RFs/{code}_MAFFT_vs_true.rfdist")
#     rf2 = extract_rf(
#         f"/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/RF_DISTANCE/{code}/RFs/{code}_MAFFT_without_SP_Col_vs_true.rfdist")
#     rf3 = extract_rf(
#         f"/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/RF_DISTANCE/{code}/RFs/{code}_indel_simulated_vs_true.rfdist")
#
#     row = {"rf1": [rf1], "rf2": [rf2], "rf3": [rf3]}
#     # row = [rf1, rf2, rf3]
#     df2 = pd.DataFrame(row)
#     df = pd.concat([df, df2], ignore_index=True)
#
#     df.to_csv("/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/RF_DISTANCE/rfs.csv", index=False)
#     print(df)

def calculate_rf():
    df = pd.DataFrame(columns=["rf1", "rf2", "rf3"])
    for fasta_path in glob.glob("/Users/kpolonsky/PycharmProjects/HoT_Py/SIMULATED_SEQ" + "/*.fas"):
        fasta_file = fasta_path.split("/")[6]
        code = fasta_file.split(".fas")[0]

        # simulated MSA
        cmd = f"cd /Users/kpolonsky/PycharmProjects/HoT_Py; rf /Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/RF_DISTANCE/{code}/Indellible_simulated_alignment/{code}_TRUE.fas.treefile /Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/RF_DISTANCE/{code}/True_tree/{code}_indelible_true_tree.ctrl"
        # subprocess.run(cmd, shell=True, check=True)
        process = subprocess.Popen(cmd, shell=True,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)

        # wait for the process to terminate
        out, err = process.communicate()
        rf1 = float(out)

        # regular MAFFT MSA
        cmd = f"cd /Users/kpolonsky/PycharmProjects/HoT_Py; rf /Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/RF_DISTANCE/{code}/MAFFT/MSA.MAFFT.aln.With_Names.treefile /Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/RF_DISTANCE/{code}/True_tree/{code}_indelible_true_tree.ctrl"
        # subprocess.run(cmd, shell=True, check=True)
        # rf2 = 0
        process = subprocess.Popen(cmd, shell=True,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)

        # wait for the process to terminate
        out, err = process.communicate()
        rf2 = float(out)

        # MSA without columns
        cmd = f"cd /Users/kpolonsky/PycharmProjects/HoT_Py; rf /Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/RF_DISTANCE/{code}/MAFFT_after_column_deletions/MSA.MAFFT.Without_low_SP_Col.With_Names.treefile /Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/RF_DISTANCE/{code}/True_tree/{code}_indelible_true_tree.ctrl"
        # subprocess.run(cmd, shell=True, check=True)
        # rf3 = 0
        process = subprocess.Popen(cmd, shell=True,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)

        # wait for the process to terminate
        out, err = process.communicate()
        rf3 = float(out)

        row = {"rf1": [rf1], "rf2": [rf2], "rf3": [rf3]}
        df2 = pd.DataFrame(row)
        df = pd.concat([df, df2], ignore_index=True)
    df.to_csv("/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/RF_DISTANCE/rfs.csv", index=False)
    print(df)


if __name__ == '__main__':
    # run_guidance_on_set(FOLDER)
    # run_msa_set_score_on_set()
    # calculate_average_auc()
    # calculate_average_auc_column_score()
    # prepare_directories_for_rf()
    calculate_rf()
    # test()
    # pass

