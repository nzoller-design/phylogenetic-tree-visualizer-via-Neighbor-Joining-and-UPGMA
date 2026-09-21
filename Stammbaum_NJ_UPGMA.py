import copy

class MethodeError(Exception):
    pass

def parse(Pfad):
    Sequenzen = {}
    aktiversequenzname = None
    aktivesequenz = ""
    with open(Pfad, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if aktiversequenzname is not None:
                    Sequenzen[aktiversequenzname] = aktivesequenz
                aktiversequenzname = line[1:].strip()
                aktivesequenz = ""
            else:
                aktivesequenz += line
        if aktiversequenzname is not None:
            Sequenzen[aktiversequenzname] = aktivesequenz
    return Sequenzen

import tkinter as tk
from tkinter import filedialog

def FastaDateiAuswählen():
    root = tk.Tk()
    root.withdraw()
    Pfad = filedialog.askopenfilename(
        title="FASTA-Datei auswählen",
        filetypes=[("FASTA-Dateien", "*.fasta *.fa *.fna"), ("Alle Dateien", "*.*")]
    )
    root.destroy()
    return Pfad

Pfad = FastaDateiAuswählen()
Seq = parse(Pfad)

def NeedlemanWunsch(Seq1, Seq2, Match=1, Mismatch=-1, Gap=-1):
    Zeilen = len(Seq1) + 1
    Spalten = len(Seq2) + 1
    Matrix = [[0] * Spalten for _ in range(Zeilen)]
    for i in range(1, Zeilen):
        Matrix[i][0] = Matrix[i - 1][0] + Gap
    for j in range(1, Spalten):
        Matrix[0][j] = Matrix[0][j - 1] + Gap
    for i in range(1, Zeilen):
        for j in range(1, Spalten):
            if Seq1[i - 1] == Seq2[j - 1]:
                Diagonal = Matrix[i - 1][j - 1] + Match
            else:
                Diagonal = Matrix[i - 1][j - 1] + Mismatch
            Oben = Matrix[i - 1][j] + Gap
            Links = Matrix[i][j - 1] + Gap
            Matrix[i][j] = max(Diagonal, Oben, Links)
    Align1 = ""
    Align2 = ""
    i, j = len(Seq1), len(Seq2)
    while i > 0 and j > 0:
        aktuellerWert = Matrix[i][j]
        if Seq1[i - 1] == Seq2[j - 1]:
            DiagonalPunkte = Match
        else:
            DiagonalPunkte = Mismatch

        if aktuellerWert == Matrix[i - 1][j - 1] + DiagonalPunkte:
            Align1 = Seq1[i - 1] + Align1
            Align2 = Seq2[j - 1] + Align2
            i -= 1
            j -= 1
        elif aktuellerWert == Matrix[i - 1][j] + Gap:
            Align1 = Seq1[i - 1] + Align1
            Align2 = "-" + Align2
            i -= 1
        else:
            Align1 = "-" + Align1
            Align2 = Seq2[j - 1] + Align2
            j -= 1

    while i > 0:
        Align1 = Seq1[i - 1] + Align1
        Align2 = "-" + Align2
        i -= 1
    while j > 0:
        Align1 = "-" + Align1
        Align2 = Seq2[j - 1] + Align2
        j -= 1

    return Align1, Align2

def VergleichsPunkte(Seq1, Seq2):
    Unterschiede = 0
    for i in range(len(Seq1)):
        if Seq1[i] != Seq2[i]:
            Unterschiede += 1
    return Unterschiede / len(Seq1)

def Distanzmatrix(Sequenzen):
    Sequenznämen = list(Sequenzen.keys())
    DMatrix = {}
    for i in range(len(Sequenznämen)):
        for j in range(i+1, len(Sequenznämen)):
            name1 = Sequenznämen[i]
            name2 = Sequenznämen[j]
            AlignSeq1, AlignSeq2 = NeedlemanWunsch(Sequenzen[name1], Sequenzen[name2])
            Punkte = VergleichsPunkte(AlignSeq1, AlignSeq2)
            if name1 not in DMatrix:
                DMatrix[name1] = {}
            if name2 not in DMatrix:
                DMatrix[name2] = {}
            DMatrix[name1][name2] = Punkte
            DMatrix[name2][name1] = Punkte
    return DMatrix
matrix = Distanzmatrix(Seq)

def NJ(DMatrix):
    DMatrix = copy.deepcopy(DMatrix)
    aktiv = set(DMatrix.keys())
    while len(aktiv) > 1:
        n = len(aktiv)
        R = {name: sum(DMatrix[name][other] for other in aktiv if other != name)
             for name in aktiv}

        besteQ = float('inf')
        Neighbours = None
        for name1 in aktiv:
            for name2 in aktiv:
                if name1 < name2:
                    Q = (n - 2) * DMatrix[name1][name2] - (R[name1] + R[name2])
                    if Q < besteQ:
                        besteQ = Q
                        Neighbours = (name1, name2)

        name1, name2 = Neighbours
        new_node = f"({name1},{name2})"
        d_ij = DMatrix[name1][name2]
        for other in aktiv:
            if other != name1 and other != name2:
                dist1 = DMatrix[name1][other]
                dist2 = DMatrix[name2][other]
                ndist = (dist1 + dist2 - d_ij) / 2
                if new_node not in DMatrix:
                    DMatrix[new_node] = {}
                DMatrix[new_node][other] = ndist
                DMatrix[other][new_node] = ndist
        aktiv.remove(name1)
        aktiv.remove(name2)
        aktiv.add(new_node)
    return list(aktiv)[0]

def UPGMA(DMatrix):
    DMatrix = copy.deepcopy(DMatrix)
    aktiv = set(DMatrix.keys())
    Groesse = {name: 1 for name in aktiv}
    while len(aktiv) > 1:
        besteDist = float('inf')
        Neighbours = None
        for name1 in aktiv:
            for name2 in aktiv:
                if name1 < name2:
                    dist = DMatrix[name1][name2]
                    if dist < besteDist:
                        besteDist = dist
                        Neighbours = (name1, name2)

        name1, name2 = Neighbours
        new_node = f"({name1},{name2})"
        groesse1 = Groesse[name1]
        groesse2 = Groesse[name2]
        for other in aktiv:
            if other != name1 and other != name2:
                dist1 = DMatrix[name1][other]
                dist2 = DMatrix[name2][other]
                ndist = (groesse1 * dist1 + groesse2 * dist2) / (groesse1 + groesse2)
                if new_node not in DMatrix:
                    DMatrix[new_node] = {}
                DMatrix[new_node][other] = ndist
                DMatrix[other][new_node] = ndist
        Groesse[new_node] = groesse1 + groesse2
        aktiv.remove(name1)
        aktiv.remove(name2)
        aktiv.add(new_node)
    return list(aktiv)[0]

def MethodeAuswählen():
    print("Welche Methode soll zur Baumkonstruktion verwendet werden?")
    print("  1: Neighbor-Joining (ungewurzelter Baum, keine Annahme einer konstanten Mutationsrate)")
    print("  2: UPGMA (gewurzelter Baum, setzt konstante Mutationsrate voraus)")
    Eingabe = input("Eingabe (1/2): ").strip()
    if Eingabe == "2":
        return "UPGMA"
    elif Eingabe == "1":
        return "NJ"
    else:
        raise MethodeError('Bitte "1" oder "2" eingeben.')

Methode = MethodeAuswählen()
if Methode == "UPGMA":
    Stammbaum = UPGMA(matrix)
else:
    Stammbaum = NJ(matrix)
print(f"Methode: {Methode}")
print(Stammbaum)

def BaumStringParsen(BaumString):
    BaumString = BaumString.strip()
    if BaumString.startswith('('):
        Tiefe = 0
        TrennIndex = None
        for i, Zeichen in enumerate(BaumString):
            if Zeichen == '(':
                Tiefe += 1
            elif Zeichen == ')':
                Tiefe -= 1
            elif Zeichen == ',' and Tiefe == 1:
                TrennIndex = i
                break
        Links = BaumString[1:TrennIndex]
        Rechts = BaumString[TrennIndex + 1:-1]
        return (BaumStringParsen(Links), BaumStringParsen(Rechts))
    else:
        return BaumString

def HöheBerechnen(Knoten):
    if isinstance(Knoten, str):
        return 0
    Links, Rechts = Knoten
    return 1 + max(HöheBerechnen(Links), HöheBerechnen(Rechts))

def BlätterZählen(Knoten):
    if isinstance(Knoten, str):
        return 1
    Links, Rechts = Knoten
    return BlätterZählen(Links) + BlätterZählen(Rechts)

def ZeichneBaum(BaumString, Dateiname="Stammbaum.png"):
    import matplotlib.pyplot as plt

    Baum = BaumStringParsen(BaumString)
    MaxHöhe = HöheBerechnen(Baum)
    AnzahlBlätter = BlätterZählen(Baum)

    fig, ax = plt.subplots(figsize=(8, max(4, AnzahlBlätter * 0.4)))
    YZähler = [0]

    def Positionieren(Knoten):
        if isinstance(Knoten, str):
            y = YZähler[0]
            YZähler[0] += 1
            x = MaxHöhe
            ax.text(x + 0.1, y, Knoten.replace("_", " "), va='center', ha='left', fontsize=10)
            return x, y
        else:
            Links, Rechts = Knoten
            x1, y1 = Positionieren(Links)
            x2, y2 = Positionieren(Rechts)
            x = MaxHöhe - HöheBerechnen(Knoten)
            y = (y1 + y2) / 2
            ax.plot([x, x], [y1, y2], color="black")
            ax.plot([x, x1], [y1, y1], color="black")
            ax.plot([x, x2], [y2, y2], color="black")
            return x, y

    Positionieren(Baum)

    ax.set_xlim(-0.5, MaxHöhe + 3)
    ax.axis("off")
    plt.tight_layout()
    plt.savefig(Dateiname, dpi=150)
    plt.show()
    print(f"Baum gespeichert als {Dateiname}")

ZeichneBaum(Stammbaum)