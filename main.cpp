//
//  main.cpp
//  RSAA

#include <iostream>
#include <iostream>
#include "la_vector.hpp"
#include <sdsl/bit_vectors.hpp>
#include "rslr2.h"

using namespace std;
using namespace sdsl;


template<typename TypeIn, typename TypeOut>
std::vector<TypeOut> read_data_binary(const std::string &filename, bool first_is_size = true, size_t max_size = std::numeric_limits<TypeIn>::max()) {
    try {
        auto openmode = std::ios::in | std::ios::binary;
        if (!first_is_size)
            openmode |= std::ios::ate;

        std::fstream in(filename, openmode);
        in.exceptions(std::ios::failbit | std::ios::badbit);

        size_t size;
        if (first_is_size)
            in.read((char *) &size, sizeof(size_t));
        else {
            size = static_cast<size_t>(in.tellg() / sizeof(TypeIn));
            in.seekg(0);
        }
        size = std::min(max_size, size);

        std::vector<TypeIn> data(size);
        in.read((char *) data.data(), size * sizeof(TypeIn));
        if constexpr (std::is_same<TypeIn, TypeOut>::value)
            return data;

        return std::vector<TypeOut>(data.begin(), data.end());
    }
    catch (std::ios_base::failure &e) {
        std::cerr << e.what() << std::endl;
        std::cerr << std::strerror(errno) << std::endl;
        exit(1);
    }
}

double createRandomInput(uint32_t setbitPercentage, uint32_t bitmapLength){
    std::random_device rd;     // only used once to initialise (seed) engine
    std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister)
    std::uniform_int_distribution<uint32_t> uniform(1,100); // guaranteed unbiased uniformly distributed integers
    inputBitmap = new bit_vector(bitmapLength,0);
    uint64_t c=0;
    for(uint32_t i=0;i<bitmapLength;i++)
        if  (uniform(rng)<=setbitPercentage)
            if ((*inputBitmap)[i]==0) {(*inputBitmap)[i]=1; c++;}
    std::ofstream out("rndBmp.bin",ios::binary);
    out.write((char*)&c,sizeof(c));
    for(uint32_t i=0;i<bitmapLength;i++)
        if ((*inputBitmap)[i]==1)
            out.write((char*)&i,sizeof(i));
    out.close();
    return (double)c/(double)bitmapLength;
}


uint64_t prepareBenchmark(char* filename){
    
    auto dataRead = read_data_binary<uint32_t, uint32_t>(filename);
    
    inputBitmap = new bit_vector(dataRead[dataRead.size()-1]+1,0);
    
    for(size_t i=0; i<dataRead.size();i++) (*inputBitmap)[dataRead[i]] = 1;
    totalsetbitCount = dataRead.size();
    
    //std::cout << "Set bit ratio :"   << (double) (dataRead.size()) / (double)(dataRead[dataRead.size()-1]+1) << std::endl;
    //std::cout << "Total bit count :" << (dataRead[dataRead.size()-1]+1) << std::endl;
    
    //std::cout <<  (double) (dataRead.size()) / (double)(dataRead[dataRead.size()-1]+1) << '\t'; //std::endl;
    //std::cout <<  (dataRead[dataRead.size()-1]+1) <<  '\t'; //std::endl;
    
    if ( 0 != inputBitmap->size() % d){
        uint64_t a = inputBitmap->size();
        inputBitmap->resize((1+(inputBitmap->size()/d))*d);
        for( ; a<inputBitmap->size();a++) (*inputBitmap)[a]=0; //pad with zeros, now the input is always a multiple of d
    }
    
    return dataRead.size();
}

void benchmark(uint64_t testsize, uint32_t bitmapLength, uint32_t bitmapSetBitsCount, char* filename, unsigned int s ){
    rrr_vector<>*                rrrInput;
    bit_vector_il<>*             ilInput;
    sd_vector<>*                 sdInput;

    rrr_vector<>::rank_1_type*   rrrRank;
    rank_support_v<1>*           v1Rank;
    rank_support_v5<>*           v5Rank;
    rank_support_il<>*           ilRank;
    sd_vector<>::rank_1_type*    sdRank;

    rrr_vector<>::select_1_type* rrrSelect;
    select_support_mcl<1>*       mclSelect;
    select_support_il<>*         ilSelect;
    sd_vector<>::select_1_type*  sdSelect;
    la_vector<uint32_t,0> *      la_opt;

    std::chrono::high_resolution_clock::time_point start;
    std::chrono::high_resolution_clock::time_point stop;
    
    cout << fixed;
    cout << setprecision(3);
    
    // prepare test queries
    uint64_t* testPosRank   = new uint64_t[testsize];
    uint64_t* testPosSelect = new uint64_t[testsize];
    std::random_device rd;     // only used once to initialise (seed) engine
    std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister)
    std::uniform_int_distribution<uint64_t> uniRank(1,bitmapLength-1); // guaranteed unbiased uniformly distributed integers between 1 and bitmapLength-1, inclusive
    std::uniform_int_distribution<uint64_t> uniSelect(1,bitmapSetBitsCount); // guaranteed unbiased uniformly distributed integers between 1 and bitmapSetBitsCount, inclusive
    for(uint64_t i=0;i <testsize;i++) {testPosRank[i] = uniRank(rng); testPosSelect[i] = uniSelect(rng);}
    
    std::cout << "R/S-scheme \t RANK in microsec. \t SELECT in microsec. \t OVERHEAD-space percentage \t" << std::endl;
    std::cout << "---------- \t ----------------- \t ------------------- \t -------------------------- \t" << std::endl;
    // test RRR
    std::cout << "RRR:" << '\t' ;
    rrrInput  = new rrr_vector<>(*inputBitmap);
    rrrRank   = new rrr_vector<>::rank_1_type(rrrInput);
    rrrSelect = new rrr_vector<>::select_1_type(rrrInput);

    uint64_t sum2=0;
    start = std::chrono::high_resolution_clock::now();
    for(uint64_t i =0; i<testsize;i++){
        sum2 += (*rrrRank)(testPosRank[i]+1);
    }
    stop = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::microseconds>( stop - start ).count() / double(testsize) << '\t' ;

    start = std::chrono::high_resolution_clock::now();
    for(uint64_t i =0; i<testsize;i++){
        sum2 += (*rrrSelect)(testPosSelect[i]+1);
    }
    stop = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::microseconds>( stop - start ).count() / double(testsize) << '\t' ;

    std::cout <<  100.00 *   (double)(size_in_bytes(*rrrInput)*8.0 + size_in_bytes(*rrrRank)*8.0 + size_in_bytes(*rrrSelect)*8.0 -  bitmapLength) /  bitmapLength << "\t" << std::endl ;
    delete rrrInput;
    delete rrrRank;
    delete rrrSelect;
    //--------------------------------------------------------
    
    
    // test SD
    std::cout << "SD:" << '\t' ;

    sdInput   = new sd_vector<>(*inputBitmap);
    sdRank    = new sd_vector<>::rank_1_type(sdInput);
    sdSelect  = new sd_vector<>::select_1_type(sdInput);

    uint64_t sum3=0;
    start = std::chrono::high_resolution_clock::now();
    for(uint64_t i =0; i<testsize;i++){
        sum3 += (*sdRank)(testPosRank[i]+1);
    }
    stop = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::microseconds>( stop - start ).count() / double(testsize) << '\t' ;

    start = std::chrono::high_resolution_clock::now();
    for(uint64_t i =0; i<testsize;i++){
        sum3 += (*sdSelect)(testPosSelect[i]+1);
    }
    stop = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::microseconds>( stop - start ).count() / double(testsize) << '\t' ;

    std::cout <<  100.00 *   (double)(size_in_bytes(*sdInput)*8.0 + size_in_bytes(*sdRank)*8.0 + size_in_bytes(*sdSelect)*8.0 -  bitmapLength) /  bitmapLength << "\t" << std::endl ;
    delete sdInput;
    delete sdRank;
    delete sdSelect;
    //--------------------------------------------------------
    
    
        
    // test V1 RANK + MCL SELECT
    std::cout << "V1+MCL:" << '\t' ;
    v1Rank    = new rank_support_v<1>(inputBitmap);
    mclSelect = new select_support_mcl<1>(inputBitmap);
    uint64_t sum5=0;
    start = std::chrono::high_resolution_clock::now();
    for(uint64_t i =0; i<testsize;i++){
        sum5 += (*v1Rank)(testPosRank[i]+1);
    }
    stop = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::microseconds>( stop - start ).count() / double(testsize) << '\t' ;

    start = std::chrono::high_resolution_clock::now();
    for(uint64_t i =0; i<testsize;i++){
        sum5 += (*mclSelect)(testPosSelect[i]+1);
    }
    stop = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::microseconds>( stop - start ).count() / double(testsize) << '\t' ;

    std::cout <<  100.00 *   (double)(size_in_bytes(*v1Rank)*8.0 + size_in_bytes(*mclSelect)*8.0 ) /  bitmapLength << "\t";
    std::cout <<  100.00 *   (double)(size_in_bytes(*v1Rank)*8.0                                 ) /  bitmapLength << "\t -> V1-space \t";
    std::cout <<  100.00 *   (double)(                             size_in_bytes(*mclSelect)*8.0 ) / bitmapLength << "\t -> MCL-space" << std::endl;
    delete v1Rank;
    //--------------------------------------------------------
    
    
    
    // test V5 RANK + MCL SELECT
    std::cout << "V5+MCL:" << '\t' ;

    v5Rank    = new rank_support_v5<>(inputBitmap);
    uint64_t sum6=0;
    start = std::chrono::high_resolution_clock::now();
    for(uint64_t i =0; i<testsize;i++){
        sum6 += (*v5Rank)(testPosRank[i]+1);
    }
    stop = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::microseconds>( stop - start ).count() / double(testsize) << '\t' ;

    start = std::chrono::high_resolution_clock::now();
    for(uint64_t i =0; i<testsize;i++){
        sum6 += (*mclSelect)(testPosSelect[i]+1);
    }
    stop = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::microseconds>( stop - start ).count() / double(testsize) << '\t' ;

    std::cout <<  100.00 *   (double)(size_in_bytes(*v5Rank)*8.0 + size_in_bytes(*mclSelect)*8.0 ) /  bitmapLength << "\t";// V5+MCL\t";
    std::cout <<  100.00 *   (double)(size_in_bytes(*v5Rank)*8.0                                 ) /  bitmapLength << " (V5-space) \t";
    std::cout <<  100.00 *   (double)(                             size_in_bytes(*mclSelect)*8.0 ) /  bitmapLength << " (MCL-space) " << std::endl;
    delete v5Rank;
    delete mclSelect;
    //--------------------------------------------------------
   
    // test INTERLEAVED
    std::cout << "IL:" << '\t' ;

    ilInput   = new bit_vector_il<>(*inputBitmap);
    ilRank    = new rank_support_il<>(ilInput);
    ilSelect  = new select_support_il<>(ilInput);

    uint64_t sum4=0;
    start = std::chrono::high_resolution_clock::now();
    for(uint64_t i =0; i<testsize;i++){
        sum4 += (*ilRank)(testPosRank[i]+1);
    }
    stop = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::microseconds>( stop - start ).count() / double(testsize) << '\t' ;

    start = std::chrono::high_resolution_clock::now();
    for(uint64_t i =0; i<testsize;i++){
        sum4 += (*ilSelect)(testPosSelect[i]+1);
    }
    stop = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::microseconds>( stop - start ).count() / double(testsize) << '\t' ;

    std::cout <<  100.00 *   (double)(size_in_bytes(*ilInput)*8.0 + size_in_bytes(*ilRank)*8.0 + size_in_bytes(*ilSelect)*8.0 -  bitmapLength) /  bitmapLength << "\t" << std::endl ;
    //--------------------------------------------------------
   
    
    unsigned int d = 256, m = 256 ;
    constructDS(d, m ,s, inputBitmap);
    std::cout << "RSAA:" << '\t' ;

    /*
    // correctness test
    for(uint64_t p=0; p<bitmapLength-1;p++){
        if (rankLR(p) !=  (*ilRank)(p+1)){
            std::cout << (*ilRank)(p+1) << std::endl;
            std::cout << rankLR(p) << std::endl;
        }
    }
    std::cout << "rank tested on all positions ..." << std::endl;
    
    for(uint64_t p=1; p<=bitmapSetBitsCount;p++){
      //  std::cout << p << std::endl;
        if (selectLR(p) !=  (*ilSelect)(p)){
            std::cout << (*ilSelect)(p) << std::endl;
            std::cout << selectLR(p) << std::endl;
        }
    }
    std::cout << "select tested on all positions ..." << std::endl;

    */
    delete ilRank;
    delete ilSelect;
    delete ilInput;
    
    //Test AA ------------------------------------------------------
    
    uint64_t sum1=0;
        start = std::chrono::high_resolution_clock::now();
        for(uint64_t i =0; i<testsize;i++){
                sum1 += rankLR(testPosRank[i]);
        }
        stop = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::microseconds>( stop - start ).count() / double(testsize) << '\t' ;
    
    
        start = std::chrono::high_resolution_clock::now();
        for(uint64_t i =0; i<testsize;i++){
            sum1 += selectLR(testPosSelect[i]);
        }
        stop = std::chrono::high_resolution_clock::now();
    
    std::cout << std::chrono::duration_cast<std::chrono::microseconds>( stop - start ).count() / double(testsize) << '\t' ;
    std::cout <<  100.00 *   (double) totalspaceusageDSinBITS/bitmapLength << "\t" << std::endl ;
  
    //---------------------------------------------------------------------
    //test LA_OPT
    std::cout << "LA(OPT):" << '\t' ;

    auto dataRead = read_data_binary<uint32_t, uint32_t>(filename);
    la_opt = new la_vector<uint32_t,0>(dataRead);

    uint64_t sum7=0;
    start = std::chrono::high_resolution_clock::now();
    for(uint64_t i =0; i<testsize;i++){
        sum7 += (*la_opt).rank(testPosRank[i]+1);
    }
    stop = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::microseconds>( stop - start ).count() / double(testsize) << '\t' ;

    start = std::chrono::high_resolution_clock::now();
    for(uint64_t i =0; i<testsize;i++){
        sum6 += (*la_opt).select(testPosSelect[i]);
    }
    stop = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::microseconds>( stop - start ).count() / double(testsize) << '\t' ;

    std::cout <<  100.00 *   (double)((*la_opt).size_in_bytes()*8.0-bitmapLength) /  bitmapLength << std::endl; 
    delete la_opt;
    //--------------------------------------------------------
   
   
}

//-----------------------------------------------------------------------------------------------------------------------------

int main(int argc, char* argv[]) {
    
    if (atoi(argv[1])==0){// test on randomly generated bitmaps
        
        
        createRandomInput(atoi(argv[2]), 1024*1024*atoi(argv[3]));

        uint64_t setbitcount = prepareBenchmark((char*)"rndBmp.bin");
        setbitRatio = (double) setbitcount / inputBitmap->size();
        
        //if ( (setbitRatio<0.2) || (setbitRatio>0.5) ) return 1;
        
        std::cout <<  "Set bit ratio: " << setbitRatio << '\t'<< std::endl;
        std::cout <<  "Input bitmap length:" << inputBitmap->size() <<  '\t' << std::endl;
        s = atoi(argv[4]);
        benchmark(1000000, inputBitmap->size(), setbitcount, (char*)"rndBmp.bin", s );

    }else{
        uint64_t setbitcount = prepareBenchmark(argv[2]);
        setbitRatio = (double) setbitcount / inputBitmap->size();
        
        //if ( (setbitRatio<0.2) || (setbitRatio>0.5) ) return 1;
        
        std::cout <<  "Set bit ratio: " << setbitRatio << '\t'<< std::endl;
        std::cout <<  "Input bitmap length:" << inputBitmap->size() <<  '\t' << std::endl;
        s = atoi(argv[3]);
        benchmark(1000000, inputBitmap->size(), setbitcount, argv[2], s );
    }

    return 0;
}
