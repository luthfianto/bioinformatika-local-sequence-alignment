from operator import itemgetter
from Bio.SubsMat.MatrixInfo import blosum62

A='MMIIPKLM'
B='AAPKLMMQ'

#fungsi untuk mendapatkan k-mers dari sequence
def kmers(string, k):
    length = len(string)                
    return [(i,string[i:i+k]) for i in range(0, length-k+1) ]

#fungsi untuk mendapatkan skor 3-mers dengan matriks blosum62 dari kedua sequence
def score(A, B):
    scores = [ blosum62.get((A[i], B[i]), blosum62.get((B[i], A[i]))) for i in range(len(A)) ]
    return sum(scores)
    
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
print(center)
centerscore=center[4]
a_index, b_index=center[0], center[1]

#iterasi dari tengah ke kiri dan ke kanan

#kiri
left=min(len( A[0:a_index+3]), len(B[0:b_index+3])); left

a_left=A[a_index-(left-3):a_index]
b_left=B[b_index-(left-3):b_index]

leftscorelist=[]
for i in range(len(a_left)):
    temp=(a_left[i], b_left[i],score(a_left[i], b_left[i]))
    if temp[2]>=0:
        leftscorelist.append(temp[2])
        left_index=i
    else:
        break
    
leftscoresum=sum(leftscorelist)
print(leftscoresum)

#kanan
right=min(len( A[a_index:]), len(B[b_index:])); right
a_right=A[a_index+3:a_index+right]
b_right=B[b_index+3:b_index+right]

rightscorelist=[]

for i in range(len(a_right)):
    temp=(a_right[i], b_right[i],score(a_right[i], b_right[i]))
    if temp[2]>=0:
        rightscorelist.append(temp[2])
        right_index=i
    else:
        break
    
sum(rightscorelist)
rightscoresum=sum(rightscorelist)
print(rightscoresum)
total_score = leftscoresum + centerscore + rightscoresum
#print(total_score)

a_kiri=a_index-len(leftscorelist)
a_kanan=a_index+3+len(rightscorelist)
b_kiri=b_index-len(leftscorelist)
b_kanan=b_index+3+len(rightscorelist)
print(A[a_kiri:a_kanan])
print(B[b_kiri:b_kanan], total_score)
#
