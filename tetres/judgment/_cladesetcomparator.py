from tetres.beast_helper.cladesetcomparator import run_cladesetcomparator
import os
import itertools
from PIL import Image
import shutil


def _cladesetcomp(cmchain, beast_applauncher, burnin=10):

    img_list = []
    cwd = os.getcwd()  # Getting current Working directory
    os.chdir(cmchain.working_dir)  # setting working dir

    if cmchain.m_MChains > 2:
        try:
            os.mkdir(f"{cmchain.working_dir}/plots/temp")
        except FileExistsError:
            raise FileExistsError("Folder temp already exists!")
        for i, j in itertools.combinations(range(cmchain.m_MChains), 2):
            run_cladesetcomparator(tree_file1=cmchain.tree_files[i],
                                   tree_file2=cmchain.tree_files[j],
                                   out_file_png=f"plots/temp/{i}_{j}_cc.png",
                                   beast_applauncher=beast_applauncher,
                                   burnin=burnin)
            img_list.append([Image.open(f"plots/temp/{i}_{j}_cc.png"), i, j])

        # print(len(img_list))
        w, h = img_list[0][0].size
        new_image = Image.new('RGB', (w * cmchain.m_MChains, h * cmchain.m_MChains),  color=(135, 135, 135))

        for img in img_list:
            new_image.paste(img[0], (img[2]*w, img[1]*h))
        new_image.save(f'plots/{cmchain.name}_clade_comp.png')

        with open('data/clade_comp.txt', 'w+') as outfile:
            for i, j in itertools.combinations(range(cmchain.m_MChains), 2):
                outfile.write(f"comparison for {i}-{j}\n\n")
                with open(f"plots/temp/{i}_{j}_cc.png.txt") as infile:
                    for line in infile:
                        outfile.write(line)
                outfile.write("-------------------------------------------------------------\n")

        shutil.rmtree(f"{cmchain.working_dir}/plots/temp")

    else:
        run_cladesetcomparator(tree_file1=cmchain.tree_files[0],
                               tree_file2=cmchain.tree_files[1],
                               out_file_png=f"plots/{cmchain.name}_clade_comp.png",
                               beast_applauncher=beast_applauncher,
                               burnin=burnin, textpath=f"data/{cmchain.name}_clade_comp.txt")

    os.chdir(cwd)  # Resetting the working directory
    return 0
