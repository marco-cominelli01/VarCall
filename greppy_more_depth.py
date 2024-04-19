#!/usr/bin/python


import sys

def filter_autosomic_recessive_variants(all_variants):
    
    # Row is a variant (it's a row in the vcf file)
    for row in all_variants:     
        
        # Splitting of the row in all its columns
        all_columns = row.split('\t')
        
        # Taking only the columns corresponding to the 3 individuals
        my_columns = all_columns[9:] 
        
        # List where to put the individuals genotypes (e.g: ["0/1","1/1","0/1"], corresponding to father, child and mother, respectively)
        fa_ch_mo = []
        
        for member in my_columns:
            if member == '.' or member in POSSIBLE_MEMBERS:
                fa_ch_mo.append('./.')   
            else:
                # In a vcf file this section is like "0/1:XXXXX", so the split is done on ':' and the element before([0]) is taken ("0/1")
                fa_ch_mo.append(member.split(":")[0])      

        father, child, mother = fa_ch_mo
        
        # SEARCH_DEPTH is explained below, in the main()
        if SEARCH_DEPTH == 'basic':         # This works as: grep "0/1.*1/1.*0/1"
            if ( (father[0]=='0' and father[2]=='1') or (father[0]=='1' and father[2]=='0') ) and (child[0]=='1' and child[2]=='1') and ( (mother[0]=='0' and mother[2]=='1') or (mother[0]=='1' and mother[2]=='0') ):
                output.append(row)
        
        elif SEARCH_DEPTH == 'standard':
            if (( (father[0]=='0' and father[2]=='1') or (father[0]=='1' and father[2]=='0') ) and (child[0]=='1' and child[2]=='1') and ( (mother[0]=='0' and mother[2]=='1') or (mother[0]=='1' and mother[2]=='0') ))
                or 
        #        output.append(row)

        else:
            print("[INVALID ARGUMENT] Please insert 'basic', 'standard' or 'high' ")

    sys.stdout.write( ('\n'.join(output)+'\n') )


def filter_autosomic_dominant_variants(all_variants):
    
    # Row is a variant (it's a row in the vcf file)
    for row in all_variants:
        
        # Splitting of the row in all its columns
        all_columns = row.split('\t')
        
        # Taking only the columns corresponding to the 3 individuals
        my_columns = all_columns[9:]
        
        # List where to put the individuals genotypes (e.g: ["0/1","1/1","0/1"], corresponding to father, child and mother, respectively)
        fa_ch_mo = []
        
        for member in my_columns:
            if member == '.' or member in POSSIBLE_MEMBERS:
                fa_ch_mo.append('./.')
            else:
                # In a vcf file this section is like "0/1:XXXXX", so the split is done on ':' and the element before([0]) is taken ("0/1")
                fa_ch_mo.append(member.split(":")[0])

        father, child, mother = fa_ch_mo

        # SEARCH_DEPTH is explained below, in the main()
        if SEARCH_DEPTH == 'basic':         # This works as: grep "0/0.*0/X.*0/0" where 'X' is any number different from 0
            if (father=='0/0') and ( (child[0]=='0' and child[2] not in ZERO_DOT) or (child[0] not in ZERO_DOT and child[2]=='0') ) and (mother=='0/0'):
                output.append(row)
        elif SEARCH_DEPTH == 'standard':
            if ((father=='0/0') and ( (child[0]=='0' and child[2] not in ZERO_DOT) or (child[0] not in ZERO_DOT and child[2]=='0') ) and (mother=='0/0')) 
                or ( father=='0/0' and (mother[0]=='0' and mother[2] not in ZERO_DOT) and ((child[0]=='0' and child[2]==mother[2]) or (child[0]==mother[2] and child[2]=='0') ))
                or mother=='0/0' and (father[0]=='0' and father[2] not in ZERO_DOT) and ((child[0]=='0' and child[2]==father[2]) or (child[0]==father[2] and child[2]=='0') )):
                output.append(row)
        else:
            print("[INVALID ARGUMENT] Please insert 'basic', 'standard' or 'high' ")

    sys.stdout.write( ('\n'.join(output)+'\n') )

    
#####################################################################################################################

input_from_vcf = sys.stdin.read()                 # Take vcf content in input

DISEASE_TYPE = sys.argv[1]                        # Either 'ar' or 'ad'

SEARCH_DEPTH = sys.argv[2]                        # Either 'standard' or 'high'

                                                    # "high" means we consider also new mutations in child (but at least one acquired from parents)
                                                    # "standard" means we consider only mutations fully inherited by parents according to the genetic pattern of disease

POSSIBLE_MEMBERS = {"father","child","mother"}
ZERO_DOT = {"0","."}

output = []

all_variants = input_from_vcf.split('\n')
all_variants = all_variants[:-1]                                # The last element is an empty string

if DISEASE_TYPE == "ar":                                        # Filtering for autosomic recessive cases
    filter_autosomic_recessive_variants(all_variants)

elif DISEASE_TYPE == 'ad':                                      # Filtering for autosomic dominant cases
    filter_autosomic_dominant_variants(all_variants)            

else:                                                           # This script handles only autosomic cases
    print("[INVALID ARGUMENT] Please insert 'ad' or 'ar' ")
    sys.exit(1)
    
