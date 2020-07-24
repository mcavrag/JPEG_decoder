// Implementation of a simple JPEG decoder
// Matija Čavrag, FER, 2016, matija.cavrag@fer.hr
// version:     2.0         2.4.2016.


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


#define EOI1 0xff
#define EOI2 0xd9
#define SOS1 0xff
#define SOS2 0xda
#define PI 3.14159


static const unsigned int q_lum[8][8] = {

{   8,   6,   5,   8,  12,  20,  26,  31 },
{   6,   6,   7,  10,  13,  29,  30,  28 },
{   7,   7,   8,  12,  20,  29,  35,  28 },
{   7,   9,  11,  15,  26,  44,  40,  31 },
{   9,  11,  19,  28,  34,  55,  52,  39 },
{  12,  18,  28,  32,  41,  52,  57,  46 },
{  25,  32,  39,  44,  52,  61,  60,  51 },
{  36,  46,  48,  49,  56,  50,  52,  50 }
};

static const unsigned int q_chrom[8][8] = {
{   9,   9,  12,  24,  50,  50,  50,  50 },
{   9,  11,  13,  33,  50,  50,  50,  50 },
{  12,  13,  28,  50,  50,  50,  50,  50 },
{  24,  33,  50,  50,  50,  50,  50,  50 },
{  50,  50,  50,  50,  50,  50,  50,  50 },
{  50,  50,  50,  50,  50,  50,  50,  50 },
{  50,  50,  50,  50,  50,  50,  50,  50 },
{  50,  50,  50,  50,  50,  50,  50,  50 },
};

FILE *fin, *fp, *fout;

char input[100*100*100*100] = {0};

double IDCTcoeff[8][8][8][8];

int dc_lum[12][2],ac_lum[16][11][2];
int dc_chrom[12][2],ac_chrom[16][11][2];

int count = 0;

int rl;

int DCvalue[3] = {0};

int dataBlock[8][8] = {0};

int idctDataBlock[8][8] = {0};

int DQFBlock[8][8] = {0};

int DRLBlock[64] = {0};

int DZZBlock[64] = {0};

float IDCTBlock[8][8] = {0.0};

int YUV[8][8][3] = {0};

int RGB[8][8][3] = {0};

int fileoffset;

char fname[128];

unsigned char outputBuf[3*512*512] = {0};

//Decoding the DC value using Huffman tables.

int getDCvalue(int * pos, int comp) {

    int flag, i, j, temp;

    int relPos, length;

    int mask, shift, huffData;

    //Searching in huffman table for DC Lum or Chrom

    for(i = 0; i < 12; i++) {

        flag = 1;

        mask = (comp == 1) ? ((1 << dc_lum[i][0])-1) : ((1 << dc_chrom[i][0]) - 1);
        shift = (comp == 1) ? (dc_lum[i][0]-1) : (dc_chrom[i][0]-1);
        huffData = (comp == 1) ? (dc_lum[i][1]) : (dc_chrom[i][1]);

        relPos = 0;

        length = (comp == 1) ? (dc_lum[i][0]):(dc_chrom[i][0]);

        for(j = length; j > 0; j--) {

            if(((huffData & mask) >> shift) != (input[*pos + relPos]-48)) {
                flag = 0;
                break;
            } else {
                relPos++;
                shift-=1;
                mask = mask >> 1;
            }
        }

        if(flag == 1) {
            break;
        }

    }

    *pos += relPos;

    temp = 0;


    for(j = 0; j < i; j++) {
        temp *=2;
        if(input[*pos+j] == '1') {
            temp++;
        }
    }

    if(input[*pos] == '0') {
        temp = temp + 1 - (int)pow(2,i);
    }

    *pos += i;

    return temp;

}


//Decoding of 63 AC values using Huffman tables.

void getACvalue(int *pos, int comp) {

    int flag, i, j, k, temp, n;

    int relPos, length;

    int mask, shift, huffData;

    int ACcounter = 1;

    while(ACcounter <= 64) {

        for(i = 0; i < 16; i++) {

            for(j = 0; j < 11; j++) {

            	flag = 1;

                mask = (comp == 1) ? ((1 << ac_lum[i][j][0])-1) : ((1 << ac_chrom[i][j][0])-1);
                shift = (comp == 1) ? (ac_lum[i][j][0]-1) : (ac_chrom[i][j][0]-1);
                huffData = (comp == 1) ? (ac_lum[i][j][1]) : (ac_chrom[i][j][1]);

                relPos = 0;

                length = (comp == 1) ? (ac_lum[i][j][0]):(ac_chrom[i][j][0]);

                for(k = length; k > 0; k--) {

                    if(((huffData & mask) >> shift) != (input[*pos + relPos]-48)) {
                        flag = 0;
                        break;
                    } else {
                        relPos++;
                        shift-=1;
                        mask = mask >> 1;
                    }   
                }


				if (input[*pos + relPos + 1] == '\0') {
					*pos += relPos;
					DRLBlock[rl++] = 0;
					DRLBlock[rl++] = 0;
					return;
				}

                if((flag == 1) && !(i != 15 && i != 0 && j == 0)) {

                	*pos += relPos;

					if ((ACcounter + i + 1) > 64) {
						DRLBlock[rl++] = 0;
						DRLBlock[rl++] = 0;
						return;
					}

                    if(i == 0 && j == 0) {
                        DRLBlock[rl++] = 0;
                        DRLBlock[rl++] = 0;

                        return;
                    }

                    n = i;
                    temp = 0;

               		if(j == 0) {
               			DRLBlock[rl++] = n;
                        DRLBlock[rl++] = temp;

               		} else {

	                    for(k = 0; k < j; k++) {
	                        temp*=2;

	                        if(input[*pos+k] == '1') {
	                            temp+=1;
	                        } 
	                    }

	                    if(input[*pos] == '0' )
	                        temp = temp + 1 - (int)pow(2, j);

	                    *pos += j;

	                    DRLBlock[rl++] = n;
	                    DRLBlock[rl++] = temp;

	                }

                    ACcounter+= i+1;

                    break;

                }

                if (ACcounter >= 64) {
                    return;
                }

                flag = 1;
            }
        }
    }
}


//The inverse run length decoding function used to retrieve AC values.

void DRLE() {
    int rld=1;
    int i=1;
    int k = 0;
    int counter = 0;
    DZZBlock[0] = DRLBlock[0];
    while(i<64) {
        if(DRLBlock[rld]==0 && DRLBlock[rld+1]==0) {
          for(k=i;k<64;k++)
            DZZBlock[k] = 0;
          return;
        }
		for (k = 0; k < DRLBlock[rld]; k++) {
			DZZBlock[i++] = 0;
		}
		if (i <= 63 && rld <= 62) {
			DZZBlock[i++] = DRLBlock[rld + 1];
			rld += 2;
		}
		else if (i == 63 && rld <= 62) {
			DZZBlock[i++] = DRLBlock[rld + 1];
		}
		else if (i <= 63 && rld == 63) {
			DZZBlock[i++] = 0;
		}
		else if (i == 63 && rld == 63) {
			DZZBlock[i] = DRLBlock[rld];
		} 
    }
}

//A simple DeZigZag function.

void DZigZag() {
    DQFBlock[0][0] = DZZBlock[0];
    DQFBlock[0][1] = DZZBlock[1];
    DQFBlock[0][2] = DZZBlock[5];
    DQFBlock[0][3] = DZZBlock[6];
    DQFBlock[0][4] = DZZBlock[14];
    DQFBlock[0][5] = DZZBlock[15];
    DQFBlock[0][6] = DZZBlock[27];
    DQFBlock[0][7] = DZZBlock[28];
    DQFBlock[1][0] = DZZBlock[2];
    DQFBlock[1][1] = DZZBlock[4];
    DQFBlock[1][2] = DZZBlock[7];
    DQFBlock[1][3] = DZZBlock[13];
    DQFBlock[1][4] = DZZBlock[16];
    DQFBlock[1][5] = DZZBlock[26];
    DQFBlock[1][6] = DZZBlock[29];
    DQFBlock[1][7] = DZZBlock[42];
    DQFBlock[2][0] = DZZBlock[3];
    DQFBlock[2][1] = DZZBlock[8];
    DQFBlock[2][2] = DZZBlock[12];
    DQFBlock[2][3] = DZZBlock[17];
    DQFBlock[2][4] = DZZBlock[25];
    DQFBlock[2][5] = DZZBlock[30];
    DQFBlock[2][6] = DZZBlock[41];
    DQFBlock[2][7] = DZZBlock[43];
    DQFBlock[3][0] = DZZBlock[9];
    DQFBlock[3][1] = DZZBlock[11];
    DQFBlock[3][2] = DZZBlock[18];
    DQFBlock[3][3] = DZZBlock[24];
    DQFBlock[3][4] = DZZBlock[31];
    DQFBlock[3][5] = DZZBlock[40];
    DQFBlock[3][6] = DZZBlock[44];
    DQFBlock[3][7] = DZZBlock[53];
    DQFBlock[4][0] = DZZBlock[10];
    DQFBlock[4][1] = DZZBlock[19];
    DQFBlock[4][2] = DZZBlock[23];
    DQFBlock[4][3] = DZZBlock[32];
    DQFBlock[4][4] = DZZBlock[39];
    DQFBlock[4][5] = DZZBlock[45];
    DQFBlock[4][6] = DZZBlock[52];
    DQFBlock[4][7] = DZZBlock[54];
    DQFBlock[5][0] = DZZBlock[20];
    DQFBlock[5][1] = DZZBlock[22];
    DQFBlock[5][2] = DZZBlock[33];
    DQFBlock[5][3] = DZZBlock[38];
    DQFBlock[5][4] = DZZBlock[46];
    DQFBlock[5][5] = DZZBlock[51];
    DQFBlock[5][6] = DZZBlock[55];
    DQFBlock[5][7] = DZZBlock[60];
    DQFBlock[6][0] = DZZBlock[21];
    DQFBlock[6][1] = DZZBlock[34];
    DQFBlock[6][2] = DZZBlock[37];
    DQFBlock[6][3] = DZZBlock[47];
    DQFBlock[6][4] = DZZBlock[50];
    DQFBlock[6][5] = DZZBlock[56];
    DQFBlock[6][6] = DZZBlock[59];
    DQFBlock[6][7] = DZZBlock[61];
    DQFBlock[7][0] = DZZBlock[35];
    DQFBlock[7][1] = DZZBlock[36];
    DQFBlock[7][2] = DZZBlock[48];
    DQFBlock[7][3] = DZZBlock[49];
    DQFBlock[7][4] = DZZBlock[57];
    DQFBlock[7][5] = DZZBlock[58];
    DQFBlock[7][6] = DZZBlock[62];
    DQFBlock[7][7] = DZZBlock[63];
}

//Dequantization of the decoded data (predefined tables for each component).

void DQuant(int comp) {
    int i, j;

    for(i = 0; i < 8; i++) {
        for(j = 0; j < 8; j++) {

            IDCTBlock[i][j] = (DQFBlock[i][j] * (float)((comp == 1) ? q_lum[i][j] : q_chrom[i][j]));
        }
    }
}

//The Inverse Discrete Cosine Transform that is used to get the YUV components

void IDCT() {
	float sum;
	int x, y, u, v;
	for (x = 0; x<8; x++)
		for (y = 0; y<8; y++) {
			sum = 0.0;
			for (v = 0; v<8; v++)
				for (u = 0; u<8; u++)
					sum += (float)(IDCTcoeff[v][u][x][y]) * IDCTBlock[v][u];
			dataBlock[x][y] = (int)(sum*(1.0 / 4.0));
		}


}

//Shift after the IDCT that is needed because of the shift that happened in DCT.

void shift() {
    int i, j;

    for(i = 0; i < 8; i++) {
        for(j = 0; j < 8; j++) {
            dataBlock[i][j] += 128;
        }
    }
}

//Storing Y,U and V components from calculated data to an array.

void store(int comp) {
    int i, j;
    for(i = 0; i < 8; i++) {
        for(j = 0; j < 8; j++) {
            YUV[i][j][comp-1] = dataBlock[i][j];
        }
    }

}

//Color conversion from YUV to RGB space.

void yuv2rgb() {
    int i, j;

   	for(i = 0; i < 8; i++) {
       	for(j = 0; j < 8; j++) {
            RGB[i][j][0] = (int)(YUV[i][j][0] + (1.402 * (YUV[i][j][2]-128)));
            RGB[i][j][1] = (int)(YUV[i][j][0] - (0.344 * (YUV[i][j][1]-128)) - (0.714 * (YUV[i][j][2]-128)));
            RGB[i][j][2] = (int)(YUV[i][j][0] + 1.772 * (YUV[i][j][1]-128));
        }
	}

    for(i = 0; i < 8; i++) {
    	for(j = 0; j < 8; j++) {
	        if(RGB[i][j][0] < 0) RGB[i][j][0] = 0;
	        else if(RGB[i][j][0] > 255) RGB[i][j][0] = 255;
	        if(RGB[i][j][1] < 0) RGB[i][j][1] = 0;
	        else if(RGB[i][j][1] > 255) RGB[i][j][1] = 255;
	        if(RGB[i][j][2] < 0) RGB[i][j][2] = 0;
	        else if(RGB[i][j][2] > 255) RGB[i][j][2] = 255;
	    }
    }

}

//Storing the RGB values of pixels in .ppm format.

void storeRGB(int x) {

	int i, j, k;

	long tmpx,tmpxx;

	tmpx=(long)x/(512/8);
	tmpxx=(long)x-tmpx*(512/8);
	fseek(fout,(long)(tmpx*(512/8)*192+tmpxx*24+fileoffset),SEEK_SET);

	for(i = 0; i < 8; i++) {
		for(j = 0; j < 8; j++) {
			for(k = 0; k < 3; k++) {
				fwrite(&RGB[i][j][k], sizeof(char), 1, fout);
			}
		}
		fseek(fout,(long)((512/8-1)*24),1);
	}	

}

//Optimization of IDCT. The cosine functions are calculated only once.

void calculateIDCTCoeff() {
	int x, y, u, v;
	for (x = 0; x<8; x++)
		for (y = 0; y<8; y++) {
			for (u = 0; u<8; u++)
				for (v = 0; v<8; v++)
					IDCTcoeff[u][v][x][y] += ((u == 0) ? 1. / sqrt(2.) : 1.) * ((v == 0) ? 1. / sqrt(2.) : 1.) *cos(((2.0*(float)x + 1.0)*(float)u*PI) / 16.0)*cos(((2.0*(float)y + 1.0)*(float)v*PI) / 16.0);
		}
}

void main(int argc, char * argv[]) {

    unsigned char data1 = 0, data2 = 0, temp = 0;

    int flag, i, j, k, position, comp, counter;

    position = 0;

	//Checking if the application was executed with right parameters.

    if (argc!=2) {
		printf("usage: decode image_file  (image_file must bi .jpg)\n");
		exit(0);
	}

	sprintf(fname, "%s.jpg", argv[1]);
	if ((fin=fopen(fname,"rb"))==NULL) {
		printf("Cannot open file: %s\n", fname);
		exit(0);
	}

	//Reading the extern Huffman tables.

    if ((fp=fopen("h_dc_lum.tbl","rb"))==NULL) {
        printf ("h_dc_lum.tbl file open error! \n");
        exit(0);
    }
    for (i=0;i<12;i++) {
        fread(&dc_lum[i][0],2,1,fp);    /// code length in bits
        fread(&dc_lum[i][1],2,1,fp); /// code word
    }
    fclose(fp);

    if ((fp=fopen("h_ac_lum.tbl","rb"))==NULL) {
        printf ("h_ac_lum.tbl file open error! \n");
        exit(0);
    }
    for (i=0;i<16;i++)                              /// i= Run
        for (j=0;j<11;j++) {                        /// j= Size
            if ((i==0) || (j>0) || (i==15)) {
                fread(&ac_lum[i][j][0],2,1,fp); /// code length
                fread(&ac_lum[i][j][1],2,1,fp); /// codeword
            }
        }
    fclose(fp);

    if ((fp=fopen("h_dc_chrom.tbl","rb"))==NULL) {
        printf ("h_dc_chrom.tbl file open error! \n");
        exit(0);
    }
    for (i=0;i<12;i++) {
        fread(&dc_chrom[i][0],2,1,fp);  /// code length in bits
        fread(&dc_chrom[i][1],2,1,fp); /// code word
    }
    fclose(fp);

    if ((fp=fopen("h_ac_chrom.tbl","rb"))==NULL) {
        printf ("h_ac_chrom.tbl file open error! \n");
        exit(0);
    }

	calculateIDCTCoeff();

    for (i=0;i<16;i++)                              // i= Run
        for (j=0;j<11;j++) {                        // j= Size
            if ((i==0) || (j>0) || (i==15)) {
                fread(&ac_chrom[i][j][0],2,1,fp);   // code length
                fread(&ac_chrom[i][j][1],2,1,fp);   // codeword
            }
        }
    fclose(fp);

    //Skipping the header for simplicity.


    while(1) {
        if(data2 == SOS1 && data1 == SOS2) {
            printf("Skipped header.\n");
            break;
        }
        data2 = data1;
        data1 = fgetc(fin);
    }

    fseek(fin, 0x0c, SEEK_CUR);

    j = 0;

	//Reading the input data from .jpg image

    while(1) {
        if(data2 == EOI1 && data1 == EOI2) {
            printf("Compressed data has been read.\n");
            break;
        }
        data2 = data1;
        data1 = fgetc(fin);
        temp = data1;

        if(data2 == 0xff && data1 == 0) {
            continue;
        }

        for(k=7 ; k >= 0 ; k--) {
            if(temp % 2 == 1) {
                input[j*8 + k] = '1';
            } else {
                input[j*8 + k] = '0';
            }

            temp/=2;
        }

        j++;
    }

    input[j*8-8] = '\0';

    comp = 1;

    counter = 0;

    flag = 1;

	//Writing header into the output file. Because of simplicity, all images are 512*512

	sprintf(fname, "%s.ppm", argv[1]);

   	if((fout = fopen(fname, "wb")) == NULL) {
		printf("Cannot open .ppm file\n");
		exit(0);
	}

    (void) fprintf(fout, "P6\n%d %d\n255\n", 512, 512);

    fileoffset=ftell(fout);

    
    while(input[position] != '\0' && flag == 1) {

        //Decoding algorithm for all input data. 
		//At the end of the loop there is another loop that checks if there are some extra bytes
		//and if there are, the decoding is done.

        rl = 0;

        DCvalue[counter] += getDCvalue(&position, comp);

        DRLBlock[rl++] = DCvalue[counter];

        getACvalue(&position, comp);

        DRLE();

        DZigZag();

        DQuant(comp);

        IDCT();

        shift();

        store(comp);

        if(counter == 2) {
            yuv2rgb();
            storeRGB(count);
            count++;
            counter = -1;
        }

        if(comp == 1) {
            comp++;
        } else if(comp == 2) {
            comp++;
        } else if (comp > 2) {
            comp=1;
        }

       	counter++;

        for(i=0;i<9;i++) {
            if(input[position+i]=='\0') {
                flag = 0;
                break;
            }
        }
    }

    fclose(fin);
   	fclose(fout);	
}