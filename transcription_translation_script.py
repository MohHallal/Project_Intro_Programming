from os.path import exists
import os
import pandas as pd

def affichage_lgn_comm(): #C'est une fonction que sera éxécuté si l'utilisateur choisit d'utiliser les lignes de commandes pour faire ses conversions
    process=input("Quel type de conversion souhaitez-vous effectuer?" "\n" "tapez 'DP' pour ADN>protéine, 'RP' pour ARN>protéine,  'DR' pour ADN>ARN:  ")
    while process != "DP" and process != "dp" and process!="RP" and process !="rp" and process!="DR" and process!="dr": # Vérification que l'utlisateur a bien suivi les consignes (si l'utilisateur utilise que les miniscules il y'aura pas de probléme)
        process = input("Erreur: Veillez à bien choisir entre 'DP', 'RP' et 'DR'!:  ")
    process2 = "N" # C'est un variable qui aura la valeur "O" si l'utilisateur veut faire la transcription et en plus la traduction, et aura la valeur "N" si l'utilisateur veut faire que la transcription (par défaut c'est "N")
    if process=="DR" or process=="dr":
        process2=input("Souhaitez-vous avoir en plus la/les séquence(s) traduite(s)? [O/N]: ")
        while process2 != "O" and process2 != "o" and process2 != "N" and process2 != "n":
            process2 = input("Erreur: Veillez à bien choisir entre 'O' pour oui et 'N' pour non! : ")
    if process=="DP" or process=="dp" or process=="RP" or process=="rp" or process2=="O" or process2=="o":  # dans ce bloc on veux savoir si l'utilisateur veux utiliser son propre code génétique
        tbl_perso=input("Souhaitez-vous utiliser une table d'a.a personalisée pour la traduction? [O/N]: ") #tbl_perso est une variable qui aura "O" si l'utilisateur veux travailler avec un code génétique personalisé, sinon il aura "N"
        while tbl_perso != "O" and tbl_perso != "o" and tbl_perso != "N" and tbl_perso != "n":
            tbl_perso = input("Erreur: Veillez à bien choisir entre 'O' et 'N': ")
        if tbl_perso=="O" or tbl_perso=="o":
            global aa_perso     # aa_perso est une variable qui contient le nom de l'acide aminé qui l'utilisateur veux lui associer un codon. ils sont mis en global pour que dans l'appel des fonctions de conversion il y'aura beacoupd d'arguments à définir
            global codon_perso  # codon_perso est une variable qui contient le codon que l'utilisateur veux l'associer a son aa_perso
            aa_perso=input("Veuillez saisir la lettre correspondant à l'a.a que vous souhaitez modifier: ")
            aa_perso=aa_perso.upper() # si l'utilisateur tape le nom de l'acide amine en miniscule, il y'aura pas de probléme
            while aa_perso=="B" or aa_perso=="J" or aa_perso=="O" or aa_perso=="U" or aa_perso=="X" or aa_perso=="Z": # Vérification de aa_perso qu'il s'agit bien d'unn acide aminé réel
                aa_perso=input("Erreur: La lettre ne correspond à aucun a.a, veillez à bien vérifier les lettres des a.a avant de saisir:  ")
            codon_perso=input("Vous souhaitez l'associer à quel codon? ")
            while not check_rna(codon_perso,"S") and len(codon_perso)!=3: # Vérification que le codon_perso est est une codon correct et contient vraiment que 3 nucléotides
                codon_perso=input("Erreur: Vous avez saisi soit un faux codon, soit un codon de plus de 3 nucléotides. Corrigez SVP! : ")
    type_seq = input("Vous dispodez de quel type de fichier?" "\n" "Tapez 'S' pour séquences, 'F' pour fichier FASTA: ") # type_seq est un variable qui aura "F" si l'utilisateur va faire les conversions à partie d'un fichier fasta,"S" si il va faire des conversions sur des séquences brutes de nucléotides
    while type_seq != "S" and type_seq != "s" and type_seq != "F" and type_seq != "f": #vérification que l'utilisateur a bien suivi les consignes
        type_seq = input("Erreur: Veillez à bien choisir entre  'S' et  'F'!: ")
    if process=="DR" or process=="dr" or process=="DP" or process=="dp": # Ce bloc sert à savoir si l'utilisateur va travailler sur une longue séquence nucléique accompagné d'un fichier GTF
        if type_seq=="F" or type_seq=="f":
            GTF=input("Est ce que votre séquence est considére comme une longue séquence nucléique (ex. un chromosome) accompagné d’un fichier GTF/GFF ? [O/N] ")
            while GTF != "O" and GTF != "o" and GTF != "N" and GTF != "n":
                GTF = input("Erreur: Veillez à bien choisir entre 'O' pour oui et 'N' pour non!  ")
            if GTF=="O" or GTF=="o":
                GTF=input("Entrez l'emplacement de votre fichier GTF: ")
                if GTF[-4:]!=".gtf" and GTF[-4:]!=".gff": # Vérification de l'extension de fichier GTF/GFF
                    print("Erreur: ce que vous avez saisi ne correspond pas à un fichier GTF.")
                    return 0
                gene=input("Entrez l'ID du gene, le nom de transcript ou du gene corespondant à la séquence que vous souhaitez convertir: ")
            else:
                del GTF # Si l'utilisateur n'a pas une longue séquence nucléique donc on aura pas besoin de ce variable, du coup on le supprime.
    if type_seq=="F" or type_seq=="f": #si l'utilisateur a un fichier FASTA, ses conversions vont étre réalisé dans ce bloc
        seq=input("Entrez l'emplacement du fichier FASTA de la séquence: ")
        if check_fasta(seq,process,type_seq): # check_fasta est une fonction qui vérifie l'existence , l'extension et le format d'un fichier FASTA, elle vérifie les nucléotides selon le type de conversion demandé, elle return True si le fichier est bon, sinon False
            if ("GTF") not in locals(): # Si l'utilisateur n'a pas une longue séquence nucléique accompagné d'un GTF, ses conversions vont étre réalisé dans ce bloc
                if process == "DR" or process == "dr": # Si l'utilisateur a choisi de faire une conversion ADN>ARN, sa conversion va étre réalisé dans ce bloc
                    result=DNA_RNA(seq,type_seq) # DNA_RNA c'est la fonction qui fait la transcription
                    RNA_result=input("Pour écupérer les résultat ( tapez 'C' pour ligne de commande, 'X' pour les exporter): ") # L'utilisateur choisit si il veux avoir ses résultats dans la ligne des commandes, ou les exporter dans un fichier qui aure le méme endroit que le fichier d'entré
                    while RNA_result != "C" and RNA_result != "c" and RNA_result != "X" and RNA_result != "x":
                        RNA_result = input("Erreur: Veillez à bien choisir entre 'C' et 'X': ")
                    if RNA_result=="C" or RNA_result=="c":
                        file_handler=open(seq,"r") # Le fichier d'entré est ouvert pour récupéré que la premiére la ligne pour que l'affichage sera plus lisible
                        if process2=="O" or process2=="o": # Si l'utilisateur a demandé la transcription est aussi la traduction, on va lui afficher sa séquence d'ARN précédé par "La séquence d'ARN:", et sa séquence d'acides aminés précédé par "La séquence protéique: "
                            print("La séquence d'ARN:","\n")
                        print(file_handler.readline(),result) # L'affichage de résultats
                        file_handler.close()
                    else: # L'exportation du résultats
                        file_handler=open(seq,"r")
                        new_directory(seq,result,file_handler.readline(),process,False) # la fonction new_directory prend l'endroit du fichier de l'entré, le résultat de la conversion, la premiére ligne à mettre en téte de résultats et crée un fichier de résultats dans le méme emplacement du fichier d'entré
                        file_handler.close()
                    if process2=="O" or process2=="o": # Si l'utilisateur a choisi de faire la transcription et la traduction , on change la valeur de processus pour que la traduction aussi réalisé
                        process="DP"
                if process=="RP" or process=="rp": # Si l'utilisateur a choisi de faire ARN>protéine, sa conversion va étre réalisé dans ce bloc
                    result = RNA_PRO(seq, type_seq,tbl_perso) # RNA_PRO est une fonction qui fait la conversion RNA>Protéine
                    PRO_result=input("Pour écupérer les résultat ( tapez 'C' pour ligne de commande, 'X' pour les exporter): ")
                    while PRO_result != "C" and PRO_result != "c" and PRO_result != "X" and PRO_result != "x":
                        PRO_result = input("Erreur: Veillez à bien choisir entre 'C' et 'X': ")
                    if PRO_result=="C" or PRO_result=="c": # Affichage de résultats dans la ligne de commandes
                        print(result)
                    else: #exportation
                        file_handler=open(seq,"r")
                        new_directory(seq,result,file_handler.readline(),process,False)
                        file_handler.close()
                if process=="DP" or process=="dp": # Si l'utilisateur a choisi de faire ADN>Protéine, sa conversion va étre utilisé dans ce bloc
                    file_handler = open(seq, "r")
                    RNA=DNA_RNA(seq,type_seq) # le variable RNA va avoir la séquence d'ARN de la séquence d'ADN d'entré, qui va étre utilisé par la fonction RNA_PRO pour produire la séquence protéique
                    directory=new_directory(seq, RNA, file_handler.readline(), "DR", True) # Cette fois on crée par la fonction new_directory un fichier FASTA provisoire qui contient la séquence d'ARN, qui va étre utilisé par RNA_PRO pour avoir le résultat attendu par l'utilisateur puis il sera supprimé (cette fonction return l'emplacement de ce fichier).
                    file_handler.close()
                    result = RNA_PRO(directory, type_seq,tbl_perso) # traduction du fichier d'ARN créé en protéine qui le résultat attendu par l'utilisateur
                    if process2=="N" or process2=="n": # Si l'utilisateur a choisi de faire la transcription et la traduction, le fichier ne sera pas supprimer.
                        os.remove(directory) # ici on supprime le fichier d'ARN créé parceque l'utilisateur a choisi d'avoir que la séquence protéique
                    PRO_result=input("Voulez vous afficher le résultat de la traduction dans la ligne de commande ou vous voulez l'exporter? taper C ou X: ")
                    while PRO_result != "C" and PRO_result != "c" and PRO_result != "X" and PRO_result != "x":
                        PRO_result = input("Erreur: Veillez à bien choisir entre 'C' et 'X'!:")
                    if PRO_result=="C" or PRO_result=="c": #Affichage des résultats
                        if process2=="O" or process2=="o":
                            print("La séquence protéique:","\n")
                        print(result)
                    else: # Exportation
                        file_handler=open(seq,"r")
                        new_directory(seq,result,file_handler.readline(),"rp",False)
                        file_handler.close()
            else: # Ce bloce va étre utilisé pour faire les conversions des longues séquences nucléique qui sont accompagnés avec des fichiers GTF
                GTF_table = pd.read_table(GTF, sep="\t", on_bad_lines='skip',header=None)
                positions=gene_positions(GTF_table,gene) # La fonction gene_positions prend le nom de géne recherché et la table de fichier GTF et fait un return de positions de start et de fin de ce géne dans le génome
                start=positions[0]
                end=positions[1]
                file_handler=open(seq,"r")
                content=file_handler.readlines()
                len_fst_line = len(content[0]) # ce variable va avoir le nombre de caractére existe dans la premiére ligne de fichier FASTA
                content="".join(content)
                start=start+len_fst_line # Le posistion de départ du géne dans le fichier
                end=end+len_fst_line # La position de fin du géne dans le fichier
                if process == "DR" or process == "dr":
                    result = DNA_RNA(content[start:end], "S")
                    print("La séquence d'ARN: ","\n" ,result)
                    if process2 == "O" or process2 == "o":
                        process = "DP"
                if process == "RP" or process == "rp":
                    result = RNA_PRO(content[start:end], "S", tbl_perso)
                    print("La séquence protéique:","\n",result)
                if process == "DP" or process == "dp":
                    RNA = DNA_RNA(content[start:end], "S")
                    result = RNA_PRO(RNA, "S", "N")
                    print("La séquence protéique: ","\n" ,result)
    else: # Ce bloc est utilisé si l'utilisateur veux faire une conversion d'une séquence brute
        seq = input("Entrez votre séquence: ")
        if process == "DR" or process == "dr":
            if check_dna(seq,type_seq): # Vérification de la séquence
                result=DNA_RNA(seq,type_seq) # Réalisation de la conversion
                print("La séquence d'ARN: ","\n",result) # affichage du résultats
            if process2=="O" or process2=="o": # Si l'utilisateur veut faire aussi la traduction en plus, on modifie le "process"
                process="rp"
        if process == "RP" or process=="rp":
            if check_rna(seq,type_seq):
                result=RNA_PRO(seq,type_seq,tbl_perso)
                print(result)
        if process=="DP" or process=="dp" or process2=="O" or process2=="o":
            if check_dna(seq,type_seq):
                RNA=DNA_RNA(seq,type_seq)
                result=RNA_PRO(RNA,type_seq,tbl_perso)
                print("La séquence protéique: ","\n",result)

def gene_positions(GTF_table,gene): #La fonction gene_positions prend le nom de géne recherché et la table de fichier GTF et fait un return de positions de start et de fin de ce géne dans le génome
    start = []
    end = []
    for i in range(1, len(GTF_table)):
        if gene in GTF_table.iloc[i, 8]:
            start.append(int(GTF_table.iloc[i, 3]))
            end.append(int(GTF_table.iloc[i, 4]))
    start=min(start)
    end=max(end)
    return start,end
def check_line_fasta(seq): # Cette fonction vérifie que tous les lignes d'un fichier FASTA ont le méme nombre de caractéres
    file_handler=open(seq,"r")
    content=file_handler.readlines()
    nbr_chr_line=content[1] # On prend le nombre de caractéres de la premiére ligne et on l'utilise pour la comparaison avec le nombre de caractéres dans le reste des lignes
    for i in range(0,len(content)-1):
        if content[i][0]!=">":
            if len(content[i])!=len(nbr_chr_line):
                if content[i+1][0]!=">" and i!=len(content)-1:
                    print("Erreur: ce fichier ne correspond pas à un FASTA, veuillez vérifier les lignes de votre fichier!")
                    file_handler.close()
                    return False
    file_handler.close()
    return True
def check_dna(seq,type_seq): # La vérification d'un fichier FASTA ou une séquence qu'il s'agit bien d'une ou plusieurs séquencs d'ADN
    if type_seq=="F" or type_seq=="f":
        ADN_list=["A","C","G","T","\n"]
        for line in seq[0:len(seq)-1]: # Ici on vériife tous les séquences du fichier FASTA sauf la derniére séquence
            for i in range(0,len(line)): # Ici on vérifie chaque caractére de chaque ligne qu'il est soit "A", "C", "G", "T", "N" (nucléotide inconnu) ou "\n" (fin de la ligne)
                if line[i] not in ADN_list:
                    print("Erreur: Séquqnce ADN fausse!")
                    return False
        drn_lgn=len(seq)-1 # ici on a le numéro d'index de la ligne ou la derniére séquence commence
        for i in seq[drn_lgn]: # Ici on vérifie la derniére séquence que c'est un ADN
            if i not in ADN_list:
                print("Erreur: Séquqnce ADN fausse!")
                print(i)
                return False
        return True
    else: # Dans ce bloc on vérifie une séquence qui est pas de type FASTA que c'est un ADN
        ADN_list = ["A", "C", "G", "T"]
        for i in seq:
            if i not in ADN_list:
                print("Erreur:séquqnce ADN fausse!")
                return False
    return True

def check_rna(seq,type_seq):
    if type_seq=="F" or type_seq=="f":
        ARN_LIST = ["A","C","G","U","\n"]
        for line in seq[0:len(seq)-1]: #dans cette boucle on va vérifie la fichier de la deuxiéme ligne jusqu'a l'avant derniére ligne, la derniére ligne n'ai pas vérifie ici parceque elle contient pas 71 caractéres
            for i in range(0,len(line)):
                if line[i] not in ARN_LIST:
                    print("Erreur: séquqnce ARN fausse")
                    return False
        drn_lgn=len(seq)-1 #dans drn_lgn j'ai mis l'index de la derniére ligne
        for i in seq[drn_lgn]: #c"est boucle for qui incrémente pour chaque caractére de la derniére ligne seulemnt, caractére par caractére
            if i not in ARN_LIST:
                print("Erreur:séquqnce ARN fausse2!")
                return False
        return True
    if type_seq == "S" or type_seq == "s":
        ARN_LIST = ["A", "C", "G", "U", "\n"]
        for i in seq:
            if i not in ARN_LIST:
                print("Erreur:séquqnce ARN fausse!")
                return False
        return True
def check_fasta(seq,process,type_seq): # check_fasta est une fonction qui vérifie l'existence , l'extension et le format d'un fichier FASTA, elle vérifie les nucléotides selon le type de conversion demandé, elle return True si le fichier est bon, sinon False
    file_existence = exists(seq)
    if file_existence: # On teste l'existence de fichier
        if seq[-3:] == ".fa" or seq[-4:] == ".fna" or seq[-4:] == ".ffn" or seq[-4:] == ".faa" or seq[-4:] == ".frn": #fa, fna, ffn, faa et frn sont tous les extensions possible d'un fichier fasta
            file_handler = open(seq, "r")
            start_fasta = file_handler.read(1)
            if start_fasta == ">": # on vérifie que le fichier commence par ">"
                if check_line_fasta(seq): # On vérifie que les lignes de fichier FASTA contient le méme nombre de caractére dans chaque ligne à part la premiére ligne
                    content=file_handler.readlines()
                    start_of_sqs=[1] # Ce variable est une liste qui contient les numéros de ligne de début des séquences. La liste est déclaré avec "1" qui est le numéro de lignede début de la premiére séquence
                    for lineno,line in enumerate(content): # Dans cette boucle on va avoir le numéro de ligne de début de chaque séquence dans un fichier FASTA
                        if line[0]==">":
                            start_of_sqs.append(lineno)
                    if process == "DR" or process == "dr" or process == "DP" or process == "dp": # Si l'utilisateur a choisi ADN>ARN ou ADN>Protéine on vérifie que sa séquence d'entré est un ADN
                        for i in range(0,len(start_of_sqs)-1):
                            if check_dna(content[start_of_sqs[i]+1:start_of_sqs[i+1]],type_seq):
                                pass
                            else:
                                file_handler.close()
                                return False
                        if check_dna(content[start_of_sqs[len(start_of_sqs)-1]+1:],type_seq):
                            pass
                        else:
                            file_handler.close()
                            return False
                        return True
                    if process == "RP" or process == "rp": # Si l'utilisateur a choisi ARN>Protéine on vérifie que la séquence d'entré est un ARN
                        for i in range(0,len(start_of_sqs)-1):
                            if check_rna(content[start_of_sqs[i]+1:start_of_sqs[i+1]],type_seq):
                                pass
                            else:
                                file_handler.close()
                                return False
                        if check_rna(content[start_of_sqs[len(start_of_sqs)-1]+1:],type_seq):
                            pass
                        else:
                            file_handler.close()
                            return False
                        return True
                else:
                    file_handler.close()
                    return False

            else:
                file_handler.close()
                print("Erreur: Cela ne correpond pas à un fichier FATSA")
                return False
        else:
            print("Erreur: D'après l'extension du fichier, cela ne correspond pas à un fichier fasta!")
            return False
    else:
        print("Erreur: Ce fichier n'existe pas!")
        return False
def DNA_RNA(seq,type_seq): # La fonction qui fait la trascription
    if type_seq=="F" or type_seq=="f": # Si l'utilisateur a un fichier FASTA
        file_handler=open(seq,"r")
        content=file_handler.readlines()
        result=""
        for line in content[1:len(content)]:
            if line[0]==">": # Si la ligne commence par > donc c'est une nouvelle séquence du coup on va l'ajouter au résultats pour faire différencier les différent séquences d'ARN a quel séquence d'ADN correspond
                result=result+line
            else:
                i = line.replace("T", "a")
                i = i.replace("A", "u")
                i = i.replace("C", "g")
                i = i.replace("G", "c")
                result=result+i.upper()
        file_handler.close()
        return result
    if type_seq=="S" or type_seq=="s": # Si l'utilisateur a une séquence d'ADN brute
        seq=seq.replace("T", "a")
        seq=seq.replace("A", "u")
        seq=seq.replace("C", "g")
        seq=seq.replace("G", "c")
        return  seq.upper()
def table_perso(table_codons):
    if codon_perso not in table_codons[aa_perso]: # Vérifie que le codon précisé par l'utilisateur n'est pas déja associé a ce acide aminé
        for aa, code in table_codons.items():
            if codon_perso in code:
                table_codons[aa].remove(codon_perso) # Suppression de l'affectation de codon sélectionné a l'acide aminé qui lui correspond naturellement
                break
        table_codons[aa_perso].append(codon_perso) # L'ajout du codon séléctionné dans les codons associé avec l'acide aminé sélectionné
    return table_codons
def RNA_PRO(seq,type_seq,tbl_perso): # La traduction
    table_codons = {'F': ['UUU', 'UUC'], 'L': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'], 'I': ['AUU', 'AUC', 'AUA'],'M': ['AUG'], 'V': ['GUU', 'GUC', 'GUA', 'GUG'], 'S': ['AGU', 'AGC', 'UCU', 'UCC', 'UCA', 'UCG'],'P': ['CCU', 'CCC', 'CCA', 'CCG'], 'T': ['ACU', 'ACC', 'ACA', 'ACG'], 'A': ['GCU', 'GCC', 'GCA', 'GCG'],'Y': ['UAU', 'UAC'], 'STOP': ['UAA', 'UAG', 'UGA'], 'H': ['CAU', 'CAC'], 'Q': ['CAA', 'CAG'],'N': ['AAU', 'AAC'], 'K': ['AAA', 'AAG'], 'D': ['GAU', 'GAC'], 'E': ['GAA', 'GAG'], 'C': ['UGU', 'UGC'],'R': ['AGA', 'AGG', 'CGU', 'CGC', 'CGA', 'CGG'], 'G': ['GGU', 'GGC', 'GGA', 'GGG']}
    if tbl_perso=="O" or tbl_perso=="o": # Changement du code génétique si l'utilisateur a demandé
        table_codons=table_perso(table_codons)
    if type_seq=="F" or type_seq=="f":
        file_handler=open(seq,"r")
        content=file_handler.readlines()
        start_of_sqs=[] # Les numéros de ligne de début de chaque séquence du fichier FASTA
        for lineno, line in enumerate(content):
            if line[0] == ">":
                start_of_sqs.append(lineno) # L'ajout des numéros de début de séquences
        result = ""
        long_ligne=len(content[1])
        for i in range(0,len(start_of_sqs)-1): # Dans ce bloc on fait la traduction de chaque séquence a part
            part=list(content[start_of_sqs[i]+1:start_of_sqs[i+1]]) # Ce variable va avoir une séquence
            result=result+content[start_of_sqs[i]]+"\n" # On ajoute la premiére ligne de la séquence pour qu'on peux identifier chaque séquence de résultat a quoi correspond
            long_fst_line=len(result) # Cette variable va prendre le nombre de caractéres présent dans chaque ligne pour que ça va étre utiliser dans le nombre de acides aminés affiché dans chaque ligne dans les résultats
            part=''.join(part)
            translate=False # Ce variable va étre True dés qu'il trouve une codon start pour initier la traduction, ou sera False dés qu'il trouve un codon stop pour arréter la traduction, par défault il est False
            for i in range(0,len(part),3):
                codon=part[i:i+3]
                if "N" in codon: # Si dans le codon il y'a un nucléotide inconnu on ajout X comme acide aminé qui ne correspond a aucun acide aminé
                    if translate: # Vérifie si la traduction est activé
                        result = result + "X"
                        if (len(result) - long_fst_line) % long_ligne == 0: # Si la ligne atteint le nombre de caractére par ligne on saute la ligne
                            result = result + "\n"
                else: # Si il n'y'a pas de nucléotides inconnu dans le codon
                    for aa, code in table_codons.items():
                        if  codon in code:
                            if aa == "M":
                                translate = True
                            if aa == "STOP":
                                translate = False
                            if translate:
                                result=result+aa
                                if (len(result)-long_fst_line)%long_ligne==0:
                                    result=result+"\n"
                                break
            result=result+"\n"
            translate=False
        part = list(content[start_of_sqs[len(start_of_sqs)-1] + 1:]) # Dans le reste de la fonction on va faire la traduction de la derniére séquence (c'est aussi cette partie seulemnt qui va étre utilisé pour faire la traduction des fichiers qui contiennent une seul séquence)
        result = result + content[start_of_sqs[len(start_of_sqs)-1]] + "\n"
        long_fst_line=len(result)
        part = ''.join(part)
        translate=False
        for i in range(0, len(part), 3):
            codon = part[i:i + 3]
            if "N" in codon:
                if translate:
                    result = result + "X"
                    if (len(result) - long_fst_line) % long_ligne == 0:
                        result = result + "\n"
            else:
                for aa, code in table_codons.items():
                    if codon in code:
                        if aa == "M":
                            translate = True
                        if aa == "STOP":
                            translate = False
                        if translate:
                            result = result + aa
                            if (len(result)-long_fst_line)%long_ligne==0:
                                result=result+"\n"
                            break
                    elif "N" in codon:
                        if translate:
                            result=result+"X"
                            if (len(result)-long_fst_line)%long_ligne==0:
                                result=result+"\n"
        return result
    else: # Le bloc qui arrive fait la traduction des séquences d'ARN brute (pas FASTA)
        result=""
        translate=False
        for i in range(0,len(seq),3):
            codon = seq[i:i + 3]
            if "N" in codon:
                if translate:
                    result = result + "X"
            else:
                for aa, code in table_codons.items():
                    if codon in code:
                        if aa == "M":
                            translate = True
                        if aa == "STOP":
                            translate = False
                        if translate:
                            result = result + aa
                            break
        return result

def new_directory(file,result,fst_line,process,temp): # la fonction new_directory prend l'endroit du fichier de l'entré, le résultat de la conversion, la premiére ligne à mettre en téte de résultats et crée un fichier de résultats dans le méme emplacement du fichier d'entré
    splices=file.split("/")
    directory=""
    for i in  range(0,len(splices)-1): # Dans ce bloc on récupére l'emplacement du fichier d'entré
        directory=directory+splices[i]+"/"
    if process=="DR" or process=="dr": # On va ajouter "RNA_of_sequence_of_" avant le nom de fichier lorsque le résultat est un ARN
        directory=directory+"RNA_of_sequence_of_"+splices[len(splices)-1]
    if process=="RP" or process=="rp": # On va ajouter "Protein_of_sequence_of_" avant le nom de fichier lorsque le résultat est une protéine
        directory = directory + "protein_of_sequence_of_" + splices[len(splices) - 1]
    new_file = open(directory, "w") # Le nouveau fichier créé
    if process!="RP" and process!="rp": # On ajoute la premiére ligne au nouveau fichier lorsque c'eest pas une conversion ARN>Protéine (parceque on a déja ajouté cette ligne au fichier pour ce type de conversion)
        new_file.write(fst_line)
    new_file.write(result) # On ajoute le résultat de conversion au fichier
    if temp==False: # Si c"est pas un fichier temporaire on affiche "fichier créé avec succée" (le fichier temporaire a été créé lorsque on veux faire la vonversion ADN>Protéine)
        print("Votre fichier a été créé avec succés!")
    new_file.close()
    return directory
if __name__ == "__main__":
    print("Bienvenue sur BIODLAD 'transcripteur/traducteur'")
    affichage_lgn_comm()
