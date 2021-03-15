#コマンドライン引数の設定
import sys
args = sys.argv

if len(args) > 1:    #引数チェック

    #args[0]にはコマンド自体が格納されているので、1から読み込む
    multi_fasta = args[1:]

    #分割されたMulti fastaファイルはfor構文でループ処理
    for i in multi_fasta:

        with open(i) as f:
            next(f)    #fastaファイルのheadを除いて読み込み
            x = f.read()

            a = x.count('a')
            t = x.count('t')
            g = x.count('g')
            c = x.count('c')

            sum = a + t + g + c
            a_p = "{:.2f}".format(100 * (a / sum))    #format関数で小数点以下の桁数を指定
            t_p = "{:.2f}".format(100 * (t / sum))
            g_p = "{:.2f}".format(100 * (g / sum))
            c_p = "{:.2f}".format(100 * (c / sum))
            gc_p = "{:.2f}".format(100 * ((g + c) / sum))

            print(i[71:] + ",A:" + str(a) + ",T:" + str(t) + ",G:" + str(g) + ",C:" + str(c) + ",A:"
                + str(a_p) + ",T:" + str(t_p) + ",G:" + str(g_p) + ",C:" + str(c_p) + ",GC:" + str(gc_p))

else:

    print("引数が指定されていません")