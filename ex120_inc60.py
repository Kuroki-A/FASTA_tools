### 60刻みで120塩基ずつ抜き出す ###

list_120 = []

with open ('File_name.fasta') as f: ### 入力File
     all_list = [l.rstrip() for l in f.readlines()]
     for i in range(len(all_list)):
         if i % 2 == 0:
             continue
         else:
            lens = len(all_list[i])
            n =  lens // 60
            if n >= 2:
                for j in range(n-1):
                    seq_tmp = all_list[i][60*j:(120+60*j)]
                    if not ('N' in  seq_tmp or 'n' in  seq_tmp):
                        list_120.append(all_list[i-1]+ '_' + str(j))
                        list_120.append(all_list[i][60*j:(120+60*j)])

with open ('File_name_cut.fasta', 'w') as g:   ### 出力File
    for i in list_120:
        print(i, file=g)
