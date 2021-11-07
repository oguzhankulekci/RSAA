//
//  rslr2.h
//  RSAA
//
#include <immintrin.h>
#include <math.h>

#ifndef rslr2_h
#define rslr2_h

typedef struct{//slope and intercept parameters of the linear regressor Y = X.slope + intercept
    float slope;
    float intercept;
}LRtype;


unsigned long long int mask64[64] = {
    0x0000000000000001, 0x0000000000000003, 0x0000000000000007, 0x000000000000000F,
    0x000000000000001F, 0x000000000000003F, 0x000000000000007F, 0x00000000000000FF,
    0x00000000000001FF, 0x00000000000003FF, 0x00000000000007FF, 0x0000000000000FFF,
    0x0000000000001FFF, 0x0000000000003FFF, 0x0000000000007FFF, 0x000000000000FFFF,
    0x000000000001FFFF, 0x000000000003FFFF, 0x000000000007FFFF, 0x00000000000FFFFF,
    0x00000000001FFFFF, 0x00000000003FFFFF, 0x00000000007FFFFF, 0x0000000000FFFFFF,
    0x0000000001FFFFFF, 0x0000000003FFFFFF, 0x0000000007FFFFFF, 0x000000000FFFFFFF,
    0x000000001FFFFFFF, 0x000000003FFFFFFF, 0x000000007FFFFFFF, 0x00000000FFFFFFFF,
    0x00000001FFFFFFFF, 0x00000003FFFFFFFF, 0x00000007FFFFFFFF, 0x0000000FFFFFFFFF,
    0x0000001FFFFFFFFF, 0x0000003FFFFFFFFF, 0x0000007FFFFFFFFF, 0x000000FFFFFFFFFF,
    0x000001FFFFFFFFFF, 0x000003FFFFFFFFFF, 0x000007FFFFFFFFFF, 0x00000FFFFFFFFFFF,
    0x00001FFFFFFFFFFF, 0x00003FFFFFFFFFFF, 0x00007FFFFFFFFFFF, 0x0000FFFFFFFFFFFF,
    0x0001FFFFFFFFFFFF, 0x0003FFFFFFFFFFFF, 0x0007FFFFFFFFFFFF, 0x000FFFFFFFFFFFFF,
    0x001FFFFFFFFFFFFF, 0x003FFFFFFFFFFFFF, 0x007FFFFFFFFFFFFF, 0x00FFFFFFFFFFFFFF,
    0x01FFFFFFFFFFFFFF, 0x03FFFFFFFFFFFFFF, 0x07FFFFFFFFFFFFFF, 0x0FFFFFFFFFFFFFFF,
    0x1FFFFFFFFFFFFFFF, 0x3FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF
};

uint64_t* rawdata;
sdsl::bit_vector* inputBitmap;
uint64_t d=256,m=256,s,logd=8,logm=8,logs=4;
double setbitRatio;
uint64_t totalsetbitCount;

sdsl::bit_vector* validityBits; // 1 bit per block to indicate whether the regression is close enough to actual value, i.e. If (0 <= (actualValue + m/2 - regressionValue) < m) them it is valid (0) else not valid(1).

unsigned char*    correctionValue; //if block is valid (validity bit is set to zero), then (corr + regressionValue - m/2) = actualValue;
LRtype*           superblockLRparams; // lineer regression parameters per each super-block
LRtype   superblockRankLR;
uint64_t totalspaceusageDSinBITS = 0, superblockCount = 0, blockCount=0;

uint64_t rankLR(uint64_t pos){
    
    uint64_t result       = 0, excess = 0;
    uint64_t offset       = pos % d;
    int64_t blockID       = pos / d;
    
    if (offset<64){
        excess = _mm_popcnt_u64(rawdata[(blockID<<2)+1]) +_mm_popcnt_u64(rawdata[(blockID<<2)+2]) + _mm_popcnt_u64(rawdata[(blockID<<2)+3]) + _mm_popcnt_u64(rawdata[(blockID<<2)] &  ~mask64[pos & 63]);
    }else if (offset<128){
        excess = (_mm_popcnt_u64(rawdata[(blockID<<2)+2]) +_mm_popcnt_u64(rawdata[(blockID<<2)+3])) +_mm_popcnt_u64(rawdata[(blockID<<2)+1] &  ~mask64[pos & 63]);
    }else if (offset<192){
        excess = _mm_popcnt_u64(rawdata[(blockID<<2)+3]) + _mm_popcnt_u64(rawdata[(blockID<<2)+2] &  ~mask64[pos & 63]);
    }else{
        excess = _mm_popcnt_u64(rawdata[(blockID<<2)+3] &  ~mask64[pos & 63]);
    }
    

    while ((blockID>=0) && ((*validityBits)[blockID]==1) && (correctionValue[blockID]==0)){
        result += _mm_popcnt_u64(rawdata[(blockID << 2)]) + _mm_popcnt_u64(rawdata[(blockID << 2)+1]) + _mm_popcnt_u64(rawdata[(blockID << 2)+2]) + _mm_popcnt_u64(rawdata[(blockID << 2)+3]);
        blockID--;
    };
    
    
    uint64_t superblockID = blockID / s;

    result += uint64_t(superblockLRparams[superblockID].intercept + superblockLRparams[superblockID].slope * (1+ blockID - superblockID*s)) - m/2 + correctionValue[blockID] - excess;
    
    if ((*validityBits)[blockID]==1){
        if (correctionValue[blockID]>=m/2)
            result += m/2;
        else
            result -= m/2;
    }

    return result;
}



uint64_t selectLR(uint64_t k){
    
    uint64_t superblockID,relativeblockID,blockID,tmp,rank=0;

    superblockID =  (uint64_t) (k/setbitRatio)/(s*d); //predict initial superblockID to start from
    while ((superblockLRparams[superblockID].intercept >= k) && (superblockID>0))  superblockID--;
    while (superblockLRparams[superblockID+1].intercept < k)  superblockID++;
    
    //predict the blockID in the predicted superblockID
    relativeblockID = (uint64_t)((k-superblockLRparams[superblockID].intercept)/superblockLRparams[superblockID].slope);
    if (relativeblockID>=s) relativeblockID = s-1;
    if (relativeblockID==0) relativeblockID = 1;
    
    blockID         = relativeblockID + (superblockID * s); //expected innerblockID
    if (blockID >= blockCount) {
        blockID =  blockCount-1;
        relativeblockID = blockID % s;
    }
    // check the total numnber of set bits till the end of the predicted block
    
    rank = (relativeblockID+1) * superblockLRparams[superblockID].slope + superblockLRparams[superblockID].intercept - m/2 + correctionValue[blockID];

    if ((*validityBits)[blockID]==1){
        if (correctionValue[blockID]==0){
            rank = rankLR(blockID*d + d-1);//detailed processing, let us call a regular rank
        }else if (correctionValue[blockID]>=m/2)
            rank += m/2;
        else
            rank -= m/2;
    }
    
    uint64_t set0 = _mm_popcnt_u64(rawdata[(blockID<<2)+0]);
    uint64_t set1 = _mm_popcnt_u64(rawdata[(blockID<<2)+1]);
    uint64_t set2 = _mm_popcnt_u64(rawdata[(blockID<<2)+2]);
    uint64_t set3 = _mm_popcnt_u64(rawdata[(blockID<<2)+3]);

    // the number of set bits in B[0..blockID*d-1] should be larger than or equal to k
    while (rank < k){
        blockID++;
        set0 = _mm_popcnt_u64(rawdata[(blockID<<2)+0]);
        set1 = _mm_popcnt_u64(rawdata[(blockID<<2)+1]);
        set2 = _mm_popcnt_u64(rawdata[(blockID<<2)+2]);
        set3 = _mm_popcnt_u64(rawdata[(blockID<<2)+3]);
        rank += set0 + set1 + set2 + set3;
    }
 
    // the number of set bits in B[0..(blockID-1)*d-1] should be smaller than  k
    rank = rank - (set0+set1+set2+set3);
    while (  rank >= k) {
        if (blockID==0)
            std::cout << "error";
        blockID--;
        set0 = _mm_popcnt_u64(rawdata[(blockID<<2)+0]);
        set1 = _mm_popcnt_u64(rawdata[(blockID<<2)+1]);
        set2 = _mm_popcnt_u64(rawdata[(blockID<<2)+2]);
        set3 = _mm_popcnt_u64(rawdata[(blockID<<2)+3]);
        rank = rank - (set0+set1+set2+set3);
    }
    
    
    tmp = k - rank ;
    
    if (set0 >= tmp )
        return (blockID << 8) + _tzcnt_u64(_pdep_u64(  1ULL << (tmp-1)   , rawdata[4*blockID] )) ;
        
    tmp -= set0 ;
    if (set1 >= tmp )
        return (blockID << 8) + 64 + _tzcnt_u64(_pdep_u64(  1ULL << (tmp-1)   , rawdata[4*blockID+1] )) ;

    
    tmp -= set1 ;
    if (set2 >= tmp )
        return (blockID << 8) + 128 + _tzcnt_u64(_pdep_u64(  1ULL << (tmp-1)   , rawdata[4*blockID+2] )) ;

    
    tmp -= set2 ;
    return (blockID << 8) + 192 + _tzcnt_u64(_pdep_u64(  1ULL << (tmp-1)   , rawdata[4*blockID+3] )) ;
                          
}
//-----------------------------------------------------------------------------------------------------------------------------


uint64_t constructDS(unsigned int d, unsigned int m, unsigned int s, sdsl::bit_vector* inputBitVector){
    
    uint64_t  setbitCount=0,errorCount=0,sel=0;
    uint64_t* selectArray=new uint64_t[(totalsetbitCount/d)+1];
    
    superblockCount = 0;
    (inputBitVector->size() % d == 0) ? blockCount = inputBitVector->size()/d : blockCount = 1 + (inputBitVector->size()/d) ;
    (blockCount % s == 0) ? superblockCount = blockCount/s : superblockCount = 1 + (blockCount/s) ;
    
    uint64_t superblockRank[superblockCount];
    superblockLRparams = new LRtype[superblockCount+1];
    
        
    uint64_t* actualValue = new uint64_t[s];
    
    validityBits      = new sdsl::bit_vector(blockCount,0);
    correctionValue   = new unsigned char[blockCount];
    
    double   SSxx=0.0,meanX = (double)(s+1)/2;
    for(unsigned int i=0; i < s; i++)  SSxx += ((double)(i+1) - meanX) * ((double)(i+1) - meanX);
    
    uint64_t sumY = 0;
    double   SSxy = 0.0, meanY;
       
    if ((*inputBitVector)[0] == 1) setbitCount=1;
    uint64_t innerBlockID = 0, superblockID = 0, blockID=0;
    
    for(uint64_t bitID=1; bitID<inputBitVector->size(); bitID++){
        
        if ((bitID % d)==0){
            
            actualValue[innerBlockID++]  = setbitCount;
            sumY += setbitCount;
            
            if ((innerBlockID % s) == 0){//compute linear regression parameters for the last s block
                
                superblockRank[superblockID] = setbitCount;
                
                meanY = (double) sumY / s;
                for(unsigned int z = 0; z < s; z++)  SSxy += ((double)(z+1) - meanX) * ((double)actualValue[z] - meanY);
                superblockLRparams[superblockID].slope     = (float) (SSxy / SSxx);
                superblockLRparams[superblockID].intercept = (float) (meanY - (meanX *  superblockLRparams[superblockID].slope));
                //test whether regression values are close enough
                for(unsigned int z = 0; z < s; z++){
                    blockID = superblockID*s+z;
                    int predictedValue = (int) (superblockLRparams[superblockID].slope * (z+1) + superblockLRparams[superblockID].intercept);
                    int residue = (int) (actualValue[z] + (m/2) - predictedValue);
                    correctionValue[blockID] = residue;

                    if ((residue<0) || (residue >= m)){
                        
                        (*validityBits)[blockID] = 1;
                        
                        if ((residue<0) && ((residue + (int)m/2)>0))
                            correctionValue[blockID] = residue + (int)m/2;
                        else if ((residue>=m) && ((residue - (int)m/2)< m))
                            correctionValue[blockID] = residue - m/2;
                        else{
                            correctionValue[blockID] = 0;
                            errorCount++;}
                    }
                }
                //---------------------------------------------------------
                superblockID++; //advance the superblockID
                innerBlockID = 0; // reset the innerblockID
                sumY         = 0; // reset the sum
                SSxy         = 0.0; // reset the SSxy
            }
        }
        if ((*inputBitVector)[bitID] == 1){
            setbitCount++;
            if (setbitCount%d==0)
                selectArray[sel++] = bitID/d;
        }
    }
    // now, fix the linear regression parameters of the last superblock which can include less than s blocks
    
    actualValue[innerBlockID++]  = setbitCount;
    sumY += setbitCount;
    meanX = (double)(innerBlockID+1)/2;
    for(unsigned int i=0; i < innerBlockID; i++)  SSxx += ((double)(i+1) - meanX) * ((double)(i+1) - meanX);
    meanY = (double) sumY / innerBlockID;
    for(unsigned int z = 0; z < innerBlockID; z++)  SSxy += ((double)(z+1) - meanX) * ((double)actualValue[z] - meanY);
    
    superblockLRparams[superblockID].slope     = (float) (SSxy / SSxx);
    superblockLRparams[superblockID].intercept = (float) (meanY - (meanX *  superblockLRparams[superblockID].slope));
    //test whether regression values are close enough
    for(unsigned int z = 0; z < innerBlockID; z++){
        blockID = superblockID*s+z;
        int predictedValue = (int) (superblockLRparams[superblockID].slope * (z+1) + superblockLRparams[superblockID].intercept);
        int residue = (int) (actualValue[z] + (m/2) - predictedValue);
        correctionValue[blockID] = residue;

        if ((residue<0) || (residue >= m)){

            (*validityBits)[blockID] = 1;
                
            if ((residue<0) && ((residue + (int)m/2)>0))
                correctionValue[blockID] = residue + (int)m/2;
            else if ((residue>=m) && ((residue - (int)m/2)< m))
                correctionValue[blockID] = residue - m/2;
            else{
                correctionValue[blockID] = 0;
                errorCount++;
            }
        }
    }
    std::cout << "(s:" << s << '\t' << "Err. Rate:" << 100.0 * (double)errorCount / (double)blockCount << ")\t" ;// std::endl;

    superblockLRparams[superblockID+1].slope = superblockLRparams[superblockID].slope;
    superblockLRparams[superblockID+1].intercept = setbitCount;
    
    rawdata = inputBitVector->data();
    
    totalspaceusageDSinBITS  = superblockCount*2*sizeof(float)*8 +
                               blockCount*8 +
                               sdsl::size_in_bytes(*validityBits)*8 +
                               2*sizeof(float)*8;
    
    return errorCount;
}




#endif /* rslr2_h */

