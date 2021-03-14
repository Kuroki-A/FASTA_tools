#やりたいこと
#１．transcript行の上にgene行をつくる
#２．gene行、transcript行の8列目先頭にそれぞれID=gene:(uniquな値), ID=transcript:(uniqueな値, geneと同値)を加える
#３．gene行8列目末尾に;biotype=protein_codingを加える
#４．transcript行8列目末尾に ID=をParent=に変えたgene行(1つ上)の8列目を加える
#５．それ以外の行(exon, CDS)の8列目先頭 Parent=gnlをParent=transcript:(uniqueな値, 属するtranscriptと同値)に変える

import pandas as pd #pandasじゃない方が簡単だったかも

df_gff3 = pd.read_table('File_name.gff3', header=None, skiprows=3) #上から3行除いて読み込み

df_transcript = df_gff3[df_gff3[2] == 'transcript'] #transcript行だけ抜き出し
df_transcript_index = df_transcript.index

#gene行のdfを作って元のdfとconcatしたいので、indexを振りなおす
new_index = [] #元のdfの新しいindex
count1 = 0
count2 = -1 #後で使う
dic = {} #後で使う

for i in range(len(df_gff3)):
    if i in df_transcript_index:
        count1 += 1
    i += count1
    new_index.append(i)

df_gff3 = pd.DataFrame(df_gff3.values, index=new_index) #元のdfのindex振り直し完了

#transcript行の8列目先頭にuniqueなIDを付け足す
df_gff3 = df_gff3.replace({8: {'ID=gnl': ''}}, regex=True) #８列目先頭削除

df_gff3.insert(3, 'add', 'ID=transcript:') #今回はIDに必要なunique値を染色体番号＋塩基配列開始位置＋終了位置にした

df_gff3[3] = '_' + df_gff3[3].astype(str) #見やすいように数字の前に_をつける
df_gff3[4] = '_' + df_gff3[4].astype(str)

df_gff3['add'] = df_gff3['add'].str.cat(df_gff3[0]) #作業用の'add'列を作りそこでIDを作る
df_gff3['add'] = df_gff3['add'].str.cat(df_gff3[3].astype(str))
df_gff3['add'] = df_gff3['add'].str.cat(df_gff3[4].astype(str))

df_gff3 = df_gff3.replace({'add': {'\.': '_'}}, regex=True) #染色体番号内にある.が気持ち悪いので_に変更
df_gff3 = df_gff3.replace({3: {'_': ''}}, regex=True) #付け足した_を消去
df_gff3 = df_gff3.replace({4: {'_': ''}}, regex=True) #付け足した_を消去

df_gff3[8][df_gff3[2]=='transcript'] = df_gff3['add'].str.cat(df_gff3[8]) #作成したIDを付け足す

df_gff3 = df_gff3.replace({'add': {'ID=transcript:': ''}}, regex=True) #transcript_id_listに渡す前処理
transcript_id_list = list(df_gff3[df_gff3[2]=='transcript']['add']) #transcript_id_listは後で使う
df_gff3 = df_gff3.drop('add', axis=1) #'add'列の削除

#gene行の作成、8列目先頭にuniqueなIDを付け足す
df_gene = df_transcript.replace({2: {'transcript': 'gene'}}) #transcript行を元にgene行を作成

gene_index =[] #indexが１つずれているので振り直し
for i in df_transcript_index:
    if i ==0:
        gene_index.append(i)
    else:
        gene_index.append(i+1)

df_gene = pd.DataFrame(df_gene.values, index=gene_index)

df_gene = df_gene.replace({8: {'ID=gnl': ''}}, regex=True) #ID付け直し

df_gene.insert(3, 'add', 'ID=gene:') #やり方はほとんどtranscript行と同じ

df_gene[3] = '_' + df_gene[3].astype(str)
df_gene[4] = '_' + df_gene[4].astype(str)

df_gene['add'] = df_gene['add'].str.cat(df_gene[0])
df_gene['add'] = df_gene['add'].str.cat(df_gene[3].astype(str))
df_gene['add'] = df_gene['add'].str.cat(df_gene[4].astype(str))

df_gene = df_gene.replace({'add': {'\.': '_'}}, regex=True)
df_gene = df_gene.replace({3: {'_': ''}}, regex=True) #付け足した_を消去
df_gene = df_gene.replace({4: {'_': ''}}, regex=True) #付け足した_を消去

df_gene[8] = df_gene['add'].str.cat(df_gene[8])
df_gene = df_gene.drop('add', axis=1) #'add'列の削除

df_gene[8] = df_gene[8] + ';biotype=protein_coding' #末尾に追加

#tanscrip行8列目の末尾にgene行の8列目を足す
df_add_transcript = df_gene.replace({8: {'ID=gene': ';Parent=gene'}}, regex=True) #IDからParentに変更、付け足し用の新しいdf作成

df_gff_add = df_gff3[8][df_gff3[2]=='transcript']
df_gff_add_index = df_gff_add.index

df_add_transcript = pd.DataFrame(df_add_transcript.values, index=df_gff_add_index) #indexをtranscript行に合わせる

df_gff3[8][df_gff3[2]=='transcript'] = df_gff3[8][df_gff3[2]=='transcript'].str.cat(df_add_transcript[8]).astype(str) #末尾に追加

#exon, CDS行のID付け直し
df_gff3 = df_gff3.replace({8: {'Parent=gnl': ''}}, regex=True)

transcript_index = df_gff3[df_gff3[2]=='transcript'].index
not_transcript_index = df_gff3[df_gff3[2]!='transcript'].index
gff3_index = df_gff3.index

for i in gff3_index:
    if i in transcript_index:
        count2 += 1
        pass
    dic[i] = transcript_id_list[count2] #exon,CDS行のindexとそれに対応するIDを辞書に保存

for i in not_transcript_index:
    df_gff3.loc[i, 8] = 'Parent=transcript:' + str(dic[i]) + '|' +  df_gff3.loc[i, 8] #Parent, ID更新

df_concat = pd.concat([df_gff3, df_gene])
df_concat = df_concat.sort_index() #concatしてindexでsort

df_concat.to_csv('File_name_format.gff3', sep='\t', header=False, index=False) #tab区切りで出力
