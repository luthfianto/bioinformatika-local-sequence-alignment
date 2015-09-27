from operator import itemgetter
from Bio.SubsMat.MatrixInfo import blosum62 #import tabel blosum62

#fungsi untuk mendapatkan k-mers dari sequence
def kmers(string, k):
    length = len(string)                
    return [(i,string[i:i+k]) for i in range(0, length-k+1) ]

#fungsi untuk mendapatkan skor 3-mers dengan matriks blosum62 dari kedua sequence
def score(A, B):
    scores = [
        #get(A,B) dan juga get(B,A), karena berbentuk tabel blosum hanya berbentuk matriks segitiga
        blosum62.get((A[i], B[i]), blosum62.get((B[i], A[i])))
        for i in range(len(A))
    ]
    
    return sum(scores)

#algoritma inti dari program blast ini    
def blastpair(A,B):
    
    #Bentuk 3-mers dari tiap sequence
    a=kmers(A,3)
    b=kmers(B,3)
    #Tentukan titik 3-mer tengah terbaik berdasar skor
    centerscorelist=[]
    for i_key,i in a:
        for j_key,j in b:
            centerscorelist.append((i_key,j_key,i,j,score(i,j)))
    
    #pilih 3-mers dengan skor tertinggi
    center=max(centerscorelist,key=itemgetter(4))
    #print(center)
    centerscore=center[4]
    a_index, b_index=center[0], center[1]
    
    ## iterasi dari tengah ke kiri dan ke kanan
    
    # iterasi kiri
    left=min(len( A[0:a_index+3]), len(B[0:b_index+3])); left
    
    a_left=A[a_index-(left-3):a_index]
    b_left=B[b_index-(left-3):b_index]
    #print(b_index)
    
    leftscorelist=[]
    for i in range(len(a_left)):
        temp=(a_left[i], b_left[i],score(a_left[i], b_left[i]))
        if temp[2]>=0:
            leftscorelist.append(temp[2])                        
        else:
            break
        
    leftscoresum=sum(leftscorelist)
    #print(leftscoresum)
    
    # iterasi kanan
    right=min(len( A[a_index:]), len(B[b_index:])); right
    a_right=A[a_index+3:a_index+right]
    b_right=B[b_index+3:b_index+right]
    
    rightscorelist=[]
    
    for i in range(len(a_right)):
        temp=(a_right[i], b_right[i],score(a_right[i], b_right[i]))
        if temp[2]>=0:
            rightscorelist.append(temp[2])
        else:
            break
        
    rightscoresum=sum(rightscorelist)
    #print(rightscoresum)
    total_score = leftscoresum + centerscore + rightscoresum
    #print(total_score)
    
    b_kiri=b_index-len(leftscorelist)
    b_kanan=b_index+3+len(rightscorelist)
    
    [print('.', end="") for i in range(b_index+2)] #print titik/spasi
    print(B[b_kiri:b_kanan], total_score, " <-", B)
    
#blastpair(A,B)

def blastdb(query, db):
    #jalankan blastpair pada tiap sequence dalam array/db
    print(query, '<- query')
    for dbsequence in db:        
        blastpair(query, dbsequence)
        
query='MMIIPKLQ' # <- contoh query
B='AAPKLMMQ' # <- contoh sequence pada db
D='AAKKLMMQ' # <- contoh sequence pada db
C=[B,D] # <- contoh db dari sequences
        
blastdb(query,C) # <- contoh memanggil blast dengan 1 query dan 1 db