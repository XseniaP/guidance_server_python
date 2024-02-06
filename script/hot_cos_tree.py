import os
import sys

from hot_cos_logger import log_print, cleanup, debug
import re
# from dendropy import Tree

class Tree_:
    def __init__(self, tree_file, user_split_number):
        self.tree_file = tree_file
        self.tree_newick = None
        self.tree_mafft = None
        self.tree_branch_split_number = user_split_number
        self.subtrees = {}
        self.splits_file = "splits.txt"

        if not os.path.isfile(f'{self.tree_file}'):
            print(f"ERROR: File not found: {self.tree_file}")
            self.tree_file = 'guide_tree.nwk'

    # Instance method
    def tree2split(self, file_path, sequence_method_command, notu, file_handler):
        br2, otu, otu2, len_val, bid = [], [], [], [], []
        otus = {}

        if not os.path.exists(file_path):
            pid = os.getpid()
            log_print(0, 2, f"ERROR: {file_handler.current_script_file} {pid} : File not found: {file_path}", file_handler)
            cleanup(1, file_handler)

        with open(file_path, 'r') as infile:
            nwstr = infile.read()

        nwstr = re.sub(r';*$', ';', nwstr)
        # nwstr = nwstr.rstrip(';').replace(' ', '')
        nwstr = nwstr.replace(' ', '')
        self.subtrees = {'tree': nwstr}

        # nwstr = nwstr[:-1]
        nwstr = nwstr.rstrip('\n')

        if notu == 2:
            if sequence_method_command == "MAF":
                match_result = re.match(r'.*:([.\d]*),.*:([.\d]*).*', nwstr)
                if match_result:
                    self.subtrees['tree'] = f"1 2 {match_result.group(1)} {match_result.group(2)}\n"
            self.subtrees['nbr'] = 0
            return self.subtrees

        nwstr = nwstr.translate(str.maketrans('()', '<>'))
        nbr = 0
        log_print(1, 0, f"nwstr=\n{nwstr}", file_handler)

        # Terminal branches
        while match_result := re.search('[<,]([^,:<>]*):([^,:<>]*)[>,]', nwstr):
            br2.append([-1, -1])
            otu.append(match_result.group(1))
            otus[match_result.group(1)] = nbr
            len_val.append(match_result.group(2))
            bid.append([nbr])
            nwstr = re.sub(f'{match_result.group(1)}:{match_result.group(2)}', str(nbr), nwstr)
            nwstr = nwstr.strip(';') + ';'
            # nwstr = nwstr[:-1]
            log_print(6, 0, f"tbrn nbr={nbr}\nnwstr=\n{nwstr}\n---\n", file_handler)
            nbr += 1

        notu = nbr

        # Internal branches
        while match_result := re.search('<(\d*),(\d*)>:([^,:<>]*)', nwstr):
            br2.append([int(match_result.group(1)), int(match_result.group(2))])
            otu.append(
                f"<{otu[int(match_result.group(1))]}:{float(len_val[int(match_result.group(1))]):.5f},{otu[int(match_result.group(2))]}:{float(len_val[int(match_result.group(2))]):.5f}>")
            len_val.append((match_result.group(3)))
            bid.append(sorted(bid[int(match_result.group(1))] + bid[int(match_result.group(2))], key=str))
            nwstr = re.sub(rf'<{match_result.group(1)},{match_result.group(2)}>:{match_result.group(3)}', str(nbr),
                           nwstr)

            log_print(6, 0, "ibrn nbr={}\nnwstr=\n{}\n---\n".format(nbr, nwstr), file_handler)
            nbr += 1

        # Last 3
        # check if rooted or unrooted
        z = re.findall(',', nwstr)
        if len(z) == 1:  # If rooted change to unrooted
            match_result = re.search(r'<(\d*),(\d*)>', nwstr)
            nbr -= 1
            i = int(match_result.group(2)) if int(match_result.group(1)) == nbr else int(match_result.group(1))
            len_val[i] = float(len_val[int(match_result.group(1))]) + float(len_val[int(match_result.group(2))])

            nwstr = nwstr.replace(str(nbr), f"{br2[nbr][0]},{br2[nbr][1]}")
            nwstr = nwstr.strip(';') + ';'
            log_print(6, 0, f"unroot nbr={nbr}\nnwstr=\n{nwstr}\n---\n", file_handler)

        otu2 = [None] * nbr
        # unrooted
        match_result = re.search('<(\d*),(\d*),(\d*)>', nwstr)
        otu2[int(match_result.group(
            1))] = f"<{otu[int(match_result.group(2))]}:{float(len_val[int(match_result.group(2))]):.5f},{otu[int(match_result.group(3))]}:{float(len_val[int(match_result.group(3))]):.5f}>"
        otu2[int(match_result.group(
            2))] = f"<{otu[int(match_result.group(3))]}:{float(len_val[int(match_result.group(3))]):.5f},{otu[int(match_result.group(1))]}:{float(len_val[int(match_result.group(1))]):.5f}>"
        otu2[int(match_result.group(
            3))] = f"<{otu[int(match_result.group(1))]}:{float(len_val[int(match_result.group(1))]):.5f},{otu[int(match_result.group(2))]}:{float(len_val[int(match_result.group(2))]):.5f}>"

        # Retrace Splits
        for i in reversed(range(nbr)):
            if br2[i][0] > -1:
                for j in range(2):
                    otu2[br2[i][j]] = f"<{otu2[i]}:{float(len_val[i]):.5f}," \
                                      f"{otu[br2[i][1 - j]]}:{float(len_val[br2[i][1 - j]]):.5f}>"

        # Recode otus
        self.subtrees['notu'] = notu
        self.subtrees['otus'] = sorted(otus.keys())

        self.subtrees['nbr'] = nbr
        self.subtrees['len'] = len_val.copy()

        i2i = [0] * notu
        oid = 0
        for oids in self.subtrees['otus']:
            i2i[int(otus[oids])] = oid
            oid += 1

        self.subtrees['br'] = []
        # Otu ids, complement, and structure fill
        splits_txt = []
        for i in range(nbr):
            a0 = [i2i[x] for x in bid[i]]
            a0 = sorted(a0)
            a1 = list(range(notu))
            for j in reversed(a0):
                a1.pop(a1.index(j))

            if len(a0) + len(a1) != notu:
                log_print(0, 2,
                          f"ERROR: subtrees error:\n  {len(a0)} : {','.join(map(str, a0))}\n  {len(a1)} : {','.join(map(str, a1))}", file_handler)
                cleanup(1, file_handler)

            self.subtrees['br'].append([{}, {}, {}])
            # subtrees['br'].append([{'notu': len(a0)}])

            self.subtrees['br'][i][0]['notu'] = len(a0)
            self.subtrees['br'][i][1]['notu'] = len(a1)
            self.subtrees['br'][i][2]['notu'] = notu
            self.subtrees['br'][i][0]['otu'] = list(a0)
            self.subtrees['br'][i][1]['otu'] = list(a1)
            self.subtrees['br'][i][2]['otu'] = list(a1) + list(a0)  # switched order on purpose, for mafft trees

            nwstr0 = otu[i].translate(str.maketrans('<>', '()'))
            self.subtrees['br'][i][0]['tree'] = nwstr0 + ";"
            nwstr = otu2[i].translate(str.maketrans('<>', '()'))
            self.subtrees['br'][i][1]['tree'] = nwstr + ";"

            if sequence_method_command == "MAF":
                otu2[int(match_result.group(
                    3))] = f"<{otu[int(match_result.group(1))]}:{float(len_val[int(match_result.group(1))]):.5f},{otu[int(match_result.group(2))]}:{float(len_val[int(match_result.group(2))]):.5f}>"
                self.subtrees['br'][i][2]['tree'] = f"({nwstr0}:" + str(
                    "{:.{}f}".format(float(len_val[i]) / 2, 5)) + f",{nwstr}:" + str(
                    "{:.{}f}".format(float(len_val[i]) / 2, 5)) + ";"
            else:
                self.subtrees['br'][i][2]['tree'] = f"({nwstr0}:{float(len_val[i]):.5f},{nwstr[1:]};"

            self.subtrees['br'][i][0]['name'] = f"b{1 if i < notu else 0}#{i:04d}"

            split_disp = f"{self.subtrees['br'][i][0]['name']} : {self.subtrees['br'][i][0]['notu']}/{self.subtrees['br'][i][1]['notu']} : [{','.join(map(str, self.subtrees['br'][i][0]['otu']))}]/[{','.join(map(str, self.subtrees['br'][i][1]['otu']))}]\n"
            splits_txt.append(split_disp)

            # debug = int(os.environ.get('DEBUG_LEVEL', 0))
            if debug > 1:
                split_disp = f"---- tree2split:\n{split_disp}  left:    {self.subtrees['br'][i][0]['tree']}\n" \
                             f"  right:    {self.subtrees['br'][i][1]['tree']}\n" \
                             f"  joined:    {self.subtrees['br'][i][2]['tree']}\n-------\n"
                log_print(6, 0, split_disp, file_handler)

        with open(self.splits_file, "w") as outfile:
            outfile.write(''.join(splits_txt))

        self.subtrees['br'].append([{}, {}, {}])
        self.subtrees['br'][int(self.subtrees['nbr'])][0]['tree'] = self.subtrees['tree']
        self.subtrees['br'][int(self.subtrees['nbr'])][0]['otu'] = list(range(self.subtrees['notu']))
        self.subtrees['br'][int(self.subtrees['nbr'])][0]['notu'] = self.subtrees['notu']
        self.subtrees['br'][int(self.subtrees['nbr'])][1]['notu'] = 0
        self.subtrees['br'][int(self.subtrees['nbr'])][2]['notu'] = 0

        self.tree_newick = self.subtrees['tree']
        # self.subtrees = self.subtrees

        if sequence_method_command == "MAF":
            self.newick2mafft(file_handler)

    def newick2mafft(self, file_handler):
        for i in range(self.subtrees['nbr'] + 1):
            for j in range(3):
                if self.subtrees['br'][i][j]['notu'] < 2:
                    continue
                log_print(6, 0, f"-br {i},{j} : {self.subtrees['br'][i][j]['tree']}\n", file_handler)
                log_print(6, 0, f"-br {i},{j} : {','.join(map(str, self.subtrees['br'][i][j]['otu']))}\n", file_handler)
                tree = self.subtrees['br'][i][j]['tree']

                # # DELETE - THIS IS TREE VISUALIZATION
                # print(tree)
                # # Create a Tree object from Newick string
                # tree__ = Tree.get_from_string(tree, schema="newick")
                # # Draw the tree
                # tree__.print_plot()

                tree = tree.translate(str.maketrans('()', '<>'))
                k = 0
                for k in range(1, self.subtrees['br'][i][j]['notu'] + 1):
                    tree = tree.replace(self.subtrees['otus'][self.subtrees['br'][i][j]['otu'][k - 1]], str(k))
                mtr = []
                while True:
                    match = re.search(r'<(\d+):([\d.]+),(\d+):([\d.]+)>', tree)
                    if match:
                        if int(match.group(1)) < int(match.group(3)):
                            k = int(match.group(1))
                            mtr.append(
                                f"{k:5d}{int(match.group(3)):5d}{float(match.group(2)):11.5f}{float(match.group(4)):11.5f}")
                        else:
                            k = int(match.group(3))
                            mtr.append(
                                f"{k:5d}{int(match.group(1)):5d}{float(match.group(4)):11.5f}{float(match.group(2)):11.5f}")
                        tree = tree.replace(match.group(0), str(k))
                    else:
                        break

                log_print(6, 0,
                          f"---\n{i} {j} : {self.subtrees['br'][i][j]['notu']} {len(mtr)}\n{self.subtrees['br'][i][j]['tree']} \n" +
                          ' ; '.join(mtr) + "\n---\n", file_handler)
                self.subtrees['br'][i][j]['tree'] = '\n'.join(mtr) + "\n"
                log_print(6, 0, f"---\n{self.subtrees['br'][i][j]['tree']}---\n", file_handler)
                # if i == 135:
                #     p=3

        log_print(3, 1, f"- {self.subtrees['tree']}\n----\n", file_handler)
        self.subtrees['tree'] = self.subtrees['br'][self.subtrees['nbr']][0]['tree']
        log_print(3, 1, f"-\n{self.subtrees['tree']}----\n", file_handler)
        self.tree_mafft = self.subtrees['tree']
        return
