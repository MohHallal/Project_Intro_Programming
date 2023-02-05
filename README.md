# Bienvenue
Ce script est un “traducteur” permettant de convertir des séquences d’ADN en séquences
d’ARN, et/ou en séquences d’acides aminés. Ce traducteur,
codé en python peux :  
– Proposer les différents niveaux de conversion : ADN > ARN, ARN > protéine, ADN >
protéine, selon le choix de l’utilisateur.  
– L’utilisateur peux fournir une ou plusieurs séquence(s) de son choix à
transcrire/traduire, ou un fichier au format FASTA contenant plusieurs séquences (le fichier FASTA est vérifié vis à vis sa: extension, nombre de caractére par ligne, chaque base).  
– L’utilisateur peux également fournir un fichier FASTA avec une longue
séquence nucléique (ex. un chromosome) accompagné d’un fichier GTF/GFF
contenant les positions des gènes d’intérêt à transcrire et/ou traduire.  
– Si le format de la séquence attendue (ex. ADN ou ARN) ne correspond pas à celle
fournie (ex. ARN), alors une erreur doit être reportée  
– Si un codon stop apparaît, alors la traduction doit s’arrêter.  
– L’utilisateur peut fournir une table d’utilisation des codons alternative des standards
(par exemple il veux Arginine pour le codon AUG et non pas le Metheonine ...).
