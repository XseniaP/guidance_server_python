import os
import shutil

from AlignmentsAddsOn import get_alignment_with_max_certainty, build_graph_from_alignments, get_tailored_alignment

if __name__ == '__main__':
    if not os.path.exists("/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/IQTREE_PYTHON_113_8_CPU_N0_convergence_0.0006/CNOT11/addsOn"):
        os.makedirs("/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/IQTREE_PYTHON_113_8_CPU_N0_convergence_0.0006/CNOT11/addsOn")
    best_alignment = get_tailored_alignment(
        path_to_folder_with_fasta_files=r"/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/IQTREE_PYTHON_113_8_CPU_N0_convergence_0.0006/CNOT11/MSA.MAFFT.Guidance2_AlternativeMSA",
        is_visualize_best_msa=True, path_to_save_png=r"/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/IQTREE_PYTHON_113_8_CPU_N0_convergence_0.0006/CNOT11/addsOn")

    print(best_alignment)
    with open(r"/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/IQTREE_PYTHON_113_8_CPU_N0_convergence_0.0006/CNOT11/addsOn/best_tailored_alignment.txt", "w") as outfile:
        for species in best_alignment['MSA']:
            outfile.write(">" + str(species)+"\n")
            outfile.write(best_alignment['MSA'][species] + "\n")

    max_msa = get_alignment_with_max_certainty(
        path_to_folder_with_fasta_files=r"/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/IQTREE_PYTHON_113_8_CPU_N0_convergence_0.0006/CNOT11/MSA.MAFFT.Guidance2_AlternativeMSA",
        is_visualize_best_msa=True, path_to_save_png=r"/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/IQTREE_PYTHON_113_8_CPU_N0_convergence_0.0006/CNOT11/addsOn")

    print(max_msa['name'])
    name = max_msa['name']
    shutil.copy(f"/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/IQTREE_PYTHON_113_8_CPU_N0_convergence_0.0006/CNOT11/MSA.MAFFT.Guidance2_AlternativeMSA/{name}", f"/Users/kpolonsky/Downloads/OrthoMaM_Simulations_GUIDANCE2/IQTREE_PYTHON_113_8_CPU_N0_convergence_0.0006/CNOT11/addsOn/{name}")
