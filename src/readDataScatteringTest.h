#include <iostream>
#include <string>
#include <fstream>
using namespace std;


// const char* fname1 = "/Users/smidt/Documents/GitHub/SCRaytrace/data/ScatteringTest/Test1.dat";
// const char* fname2 = "/Users/smidt/Documents/GitHub/SCRaytrace/data/ScatteringTest/Test2.dat";


//Converts an unsigned char representation of a byte array to an integer
int reverseIntBytes(unsigned char* bytes){
    unsigned char* bytes2;
    #if BIG_ENDIAN
        bytes2 = new unsigned char[4];
        // *fsize = (bytes[0] << 24) + (bytes[1] << 16) + (bytes[2] << 8) + bytes[3];
        bytes2[0] = bytes[3];
        bytes2[1] = bytes[2];
        bytes2[2] = bytes[1];
        bytes2[3] = bytes[0];
    #else
        // *fsize = bytes[0] + (bytes[1] << 8) + (bytes[2] << 16) + (bytes[3] << 24);
        bytes2 = bytes;
    #endif
    int result = *(int*)bytes2;
    #if BIG_ENDIAN
        delete[] bytes2;
    #endif

    return result;
}


//Converts an unsigned char representation of a byte array to a float
float reverseFloatBytes(unsigned char* bytes){
    unsigned char* bytes2;
    #if BIG_ENDIAN
        bytes2 = new unsigned char[4];
        // *fsize = (bytes[0] << 24) + (bytes[1] << 16) + (bytes[2] << 8) + bytes[3];
        bytes2[0] = bytes[3];
        bytes2[1] = bytes[2];
        bytes2[2] = bytes[1];
        bytes2[3] = bytes[0];
    #else
        // *fsize = bytes[0] + (bytes[1] << 8) + (bytes[2] << 16) + (bytes[3] << 24);
        bytes2 = bytes;
    #endif
    float result = *(float*)bytes2;
    #if BIG_ENDIAN
        delete[] bytes2;
    #endif

    return result;
}


//Reads and parses data from a .dat file and returns a float pointer
float* readDatFile(const char* fname){
    ifstream file;
    file.open(fname, ios_base::in | ios_base::binary);

    int fsize;
    unsigned char bytes[4]; // Array to store bytes

    file.read(reinterpret_cast<char*>(&bytes), sizeof(int));
    
    fsize = reverseIntBytes(bytes);

    float* realData = (float*)calloc(fsize*fsize, sizeof(float));

    for(int i = 0; i < fsize*fsize; i++){
        file.read(reinterpret_cast<char*>(&bytes), sizeof(float));
        realData[i] = reverseFloatBytes(bytes);
    }

    // printf("%d\n", fsize);
    // for(int i = 0; i < fsize*fsize; i++){
    //     printf("%.20f\n", realData[i]);
    // }

    file.close();
    return realData;
}


// int main(){
//     // float* data1 = readDatFile(fname1);
//     // float* data2 = readDatFile(fname2);

//     // free(data1);
//     // free(data2);

//     return EXIT_SUCCESS;

// }