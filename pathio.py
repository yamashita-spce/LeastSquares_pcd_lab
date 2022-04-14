#python IOにexplorerを使用する
#2021/06

import tkinter as tk
import tkinter.filedialog as fd

#explorerからファイル、フォルダーパスを取得

#インポートファイルのpathの取得
def infle():
    root = tk.Tk()
    root.withdraw()
    inf = fd.askopenfilename()
    return(inf)

#アウトプットファイルを選択・作成し、そのpathを取得
def outfle():
    root = tk.Tk()
    root.withdraw()
    ouf = fd.asksaveasfilename()
    return(ouf)

#ディレクトリのパスを取得
def indir():
    root = tk.Tk()
    root.withdraw()
    ind = fd.askdirectory()
    return(ind)
    
