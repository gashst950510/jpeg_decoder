#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define UC unsigned char 
#define PI 3.14159265358979323846

int zigzag[8][8] = { //zigzag order
	{ 0,  1,  5,  6, 14, 15, 27, 28},
	{ 2,  4,  7, 13, 16, 27, 29, 42},
	{ 3,  8, 12, 17, 25, 30, 41, 43},
	{ 9, 11, 18, 24, 31, 40, 44, 53},
	{10, 19, 23, 32, 39, 45, 52, 54},
	{20, 22, 33, 38, 46, 51, 55, 60},
	{21, 34, 37, 47, 50, 56, 59, 61},
	{35, 36, 48, 49, 57, 58, 62, 63}
};

typedef struct s_Node{ //Node for Huffman tree
	struct s_Node* left;
	struct s_Node* right;
	int value;
}Node;

struct s_JPEG{ //struct that store all needed data
	int X, Y; //width & height of picture
	
	int Nf; //number of frame, assume 3
	int Hmax, Vmax; //max on sample factor
	int H[3], V[3]; //sample factor
	int Tdi[3], Tai[3], Tqi[3]; //table index
	int DC[3]; //DC value for DPCM
	
	double cosine[8][8]; //cosine for IDCT
	double factor[8][8]; //factor for IDCT
	
	int DQT[4][64]; //id, 8x8 array
	Node* DHT[4][2]; //id, 0:DC/1:AC
	
	int** YCbCr[3]; //YCbCr data
	int** RGB[3]; //RGB data
}JPEG;

Node* new_node(){ //initialize a new node
	Node* root = malloc( sizeof(Node) );
	root->left = NULL;
	root->right = NULL;
	root->value = -1; //init to -1
	return root;
}

Node* insert_node(Node* root, int code, int len, int value){ //insert value into tree
	Node* cur = root;
	for(int mask = 1<<(len-1); mask!=0; mask>>=1){
		if(code&mask){ //left
			if(cur->left == NULL) //not yet passed, initialize
				cur->left = new_node();
			cur = cur->left;
		}
		else{ //right
			if(cur->right == NULL)
				cur->right = new_node();
			cur = cur->right;
		}
	}
	cur->value = value;
	return root;
}

int read_byte(FILE* fp){ //read a byte from jpeg, assume no Rst
	UC buffer;
	fread( &buffer, sizeof(UC), 1, fp);

	if(buffer == (UC)0xFF){ //exception when get 0xFF
		while(buffer == (UC)0xFF) //keep reading 0xFF if next byte is still 0xFF
			fread( &buffer, sizeof(UC), 1, fp);
		
		if(buffer == (UC)0x00) //if 0xFF00, return 0xFF
			return 0xFF;
		else //else, return non-0xFF code
			return buffer;		
	}
	return buffer;
}

void prepare_IDCT(){ //preparation before IDCT
	for(int n1=0; n1<8; n1++)
		for(int k1=0; k1<8; k1++)
			JPEG.cosine[n1][k1] = cos((double)((2*n1+1) * k1 * PI / 16));

	for(int k1=0; k1<8; k1++){
		for(int k2=0; k2<8; k2++){
			double c1 = k1==0 ? 1/sqrt(2):1;
			double c2 = k2==0 ? 1/sqrt(2):1;
			JPEG.factor[k1][k2] = c1*c2/4;
		}
	}
	return;
}

int get_DC(int id, FILE* fp, int* code, int* counter){ //get DC value from Huffman tree
	int msk = 1<<7; //only get the first bit of one byte
	Node* DC = JPEG.DHT[id][0]; //root of DC Huffman tree
	while(DC->value == -1){ //decode until find a value not equal to -1
		if(msk & *code) DC = DC->left;
		else DC = DC->right;
		*code <<= 1;
		(*counter)++;
		if(*counter == 8){ //when counter equals to 8, get a new byte
			*code = read_byte(fp);
			*counter = 0;
		}
	}
	
	int get_len = DC->value;
	int diff = 0;
	int sign = (msk & *code) ? 1:-1; //-1 if first bit is 0
	for(int i=0; i<get_len; i++){
		if(msk & *code) diff+=1;
		*code <<= 1;
		(*counter)++;
		if((*counter) == 8){
			*code = read_byte(fp);
			*counter = 0;
		}
		if(i != get_len-1) diff <<= 1;
	}
	
	if(sign == -1)
		diff = ((~diff) & ((1<<get_len)-1) ) * sign; //if first bit is 0, reslut will be bitwise complement times -1
	 
	return diff;
}

int get_AC(int id, FILE* fp, int* code, int* counter, int** array, int* num){ //get AC value from Huffman tree
	int msk = 1<<7; //only get the first bit of one byte
	Node* AC = JPEG.DHT[id][1]; //root of AC Huffman tree
	while(AC->value == -1){ //decode until find a value not equal to -1
		if(msk & *code) AC = AC->left; 
		else AC = AC->right; 
		*code <<= 1;
		(*counter)++;
		if(*counter == 8){
			*code = read_byte(fp);
			*counter = 0;
		}
	}
	
	int num_zero = (AC->value & 0xF0)>>4; //first 4 bit of decoded value
	for(int i=0; i<num_zero; i++){ //add some zeros
		**array = 0;
		(*array)++;
		(*num)++;
	}
	
	int get_len = AC->value & 0x0F; //last 4 bit of decoded value
	int diff = 0;
	int sign = (msk & *code) ? 1:-1;
	for(int i=0; i<get_len; i++){
		if(msk & *code) diff+=1;
		*code <<= 1;
		
		(*counter)++; 
		if(*counter == 8){
			*code = read_byte(fp);
			*counter = 0;
		}
		if(i != get_len-1) diff <<= 1;
	}
	if(sign == -1)
		diff = ((~diff) & ((1<<get_len)-1) ) * sign;
	
	**array = diff; //put the value to last element
	return AC->value;
}

int get_block(int array[3][16][16], FILE* fp, int* code, int* counter){ //get a basic block
	int buffer[64];
	int buffer1[64];
	double buffer2[8][8];
	for(int c=0; c<3; c++){
		int Tdi = JPEG.Tdi[c];
		int Tai = JPEG.Tai[c];
		int Tqi = JPEG.Tqi[c];
		
		int ratioY = JPEG.Vmax / JPEG.V[c]; //how many rows form a basic block
		int ratioX = JPEG.Hmax / JPEG.H[c]; //how many columns form a basic block

		for(int i=0; i<JPEG.V[c]; i++){
			for(int j=0; j<JPEG.H[c]; j++){
				memset( buffer, 0, 8*8 * sizeof(int)); //set 0 to a 8*8
				memset( buffer2, 0, 8*8*sizeof(double));
				buffer[0] = get_DC(Tdi, fp, code, counter) + JPEG.DC[c]; //DPCM
				JPEG.DC[c] = buffer[0];

				int num = 1; //number of elements already read in to buffer
				int* ptr = buffer+1; 
				while(get_AC(Tai, fp, code, counter, &ptr, &num) != 0){
					ptr++;
					num++;
					if(num >= 64) break; //if already read 64 elements, break
				}			
				
				//apply DQT
				int* dqt =JPEG.DQT[Tqi];
				for(int i=0; i<64; i++)
					buffer[i] *= dqt[i];
				
				//apply zigzag
				for(int i=0; i<8; i++){
					for(int j=0; j<8; j++){
						int id = zigzag[i][j];
						int pos = i*8 + j;
						buffer1[pos] = buffer[id];
					}
				}

				//apply IDCT
				for(int n1=0; n1<8; n1++)
					for(int n2=0; n2<8; n2++)
						for(int k1=0; k1<8; k1++)
							for(int k2=0; k2<8; k2++)
								buffer2[n1][n2] += JPEG.factor[k1][k2] * JPEG.cosine[n1][k1] * JPEG.cosine[n2][k2] * (double)buffer1[k1*8 + k2];	
							
				//combine small array, shift 128
				int offsetY = i * 8;
				int offsetX = j * 8;
				for(int y=0; y<8; y++)
					for(int x=0; x<8; x++)
						array[c][offsetY + y*ratioY][offsetX + x*ratioX] = (int)(round(buffer2[y][x]) + 128); 
			}
		}
		//shift & copy (recover from subsampling)
		if(JPEG.V[c] != JPEG.Vmax){
			for(int y=0; y<8; y++){
				int row = y*ratioY + 1;
				for(int x=0; x<8*ratioX; x++)
					array[c][row][x] = array[c][row-1][x];
			}
		}
		if(JPEG.H[c] != JPEG.Hmax){
			for(int x=0; x<8; x++){
				int col = x*ratioX + 1;
				for(int y=0; y<8*ratioY; y++)
					array[c][y][col] = array[c][y][col-1];
			}
		}
	}
}

void scan(FILE* fp){ //scan data code
	int code = read_byte(fp);
	int counter = 0; //count how many bit is shifted
	int Y = JPEG.Y;
	int X = JPEG.X; 
	int perY = JPEG.Vmax * 8; //height of a basic block
	int perX = JPEG.Hmax * 8; //width of a basic block
	int nY = ceil(Y / perY); //number of basic blocks in row
	int nX = ceil(X / perX); //number of basi blocks in column
	if(Y % perY != 0) nY += 1; //add one if not enough
	if(X % perX != 0) nX += 1;

	for(int i=0; i<3; i++){ 
		JPEG.DC[i] = 0; //DC value for DPCM
		//malloc space to store data information
		JPEG.YCbCr[i] = malloc( nY * perY * sizeof(int*)); 
		for(int j=0; j<(nY+1) * perY; j++)
			JPEG.YCbCr[i][j] = malloc( nX * perX * sizeof(int));
	}
	int array[3][16][16] = {0}; //array to store one basic block
	
	prepare_IDCT(); //prepare factor and cosine for IDCT
	
	for(int i=0; i<nY; i++){
		for(int j=0; j<nX; j++){
			get_block(array, fp, &code, &counter); //get a block, store in array

			//merge all into one
			int offsetY = i * perY; //offset when placing into array
			int offsetX = j * perX;
			for(int c=0; c<3; c++)
				for(int y=0; y<perY; y++)
					for(int x=0; x<perX; x++)
						JPEG.YCbCr[c][y + offsetY][x + offsetX] = array[c][y][x]; 
		}
	}
	return;
}

void YCbCr_to_RGB(){ //compute RGB data based on YCbCr
	int Y = JPEG.Y;
	int X = JPEG.X;
	//malloc size equal to picture size
	for(int i=0; i<3; i++){ 
		JPEG.RGB[i] = malloc(Y * sizeof(int*));
		for(int j=0; j<Y; j++)
			JPEG.RGB[i][j] = malloc(X * sizeof(int));
	}

	for(int i=0; i<Y; i++){
		for(int j=0; j<X; j++){
			//value of Y, Cb, Cr
			int Y = JPEG.YCbCr[0][i][j];
			int Cb = JPEG.YCbCr[1][i][j];
			int Cr = JPEG.YCbCr[2][i][j];
			//translate using equation
			JPEG.RGB[0][i][j] = round(1.402*(Cr-128) + Y);
			JPEG.RGB[1][i][j] = round(-0.34414*(Cb-128)-0.71414*(Cr-128) + Y);
			JPEG.RGB[2][i][j] = round(1.772*(Cb-128) + Y);
			//limit range (0~255)
			for(int c=0; c<3; c++){
				if(JPEG.RGB[c][i][j] < 0) 
					JPEG.RGB[c][i][j] = 0;
				if(JPEG.RGB[c][i][j] > 255) 
					JPEG.RGB[c][i][j] = 255;
			}
		}
	}
	return;
}

void write_bmp(char* filename){ //write RGB data to .bmp file, modified based on https://stackoverflow.com/questions/2654480/writing-bmp-image-in-pure-c-c-without-other-libraries
	char outfile[50];
	int k;
	for(k=0; filename[k]!='.'; k++) outfile[k] = filename[k]; //make output filename
	strcpy(outfile+k, ".bmp\0");
	
	int w = JPEG.X; //width of picture
	int h = JPEG.Y; //height of picture
	int filesize = 54 + 3*w*h; 
	UC* img  = malloc(3*w*h);
	memset(img, 0, 3*w*h);
	
	//data of bmp
	for(int i=0; i<w; i++){
		for(int j=0; j<h; j++){
			int id = (i + j*w) * 3;
			img[id+2] = (UC)JPEG.RGB[0][j][i];
			img[id+1] = (UC)JPEG.RGB[1][j][i];
			img[id+0] = (UC)JPEG.RGB[2][j][i];
		}
	}
	//header of bmp
	UC bmpfileheader[14] = {'B','M', 0, 0, 0, 0, 0, 0, 0, 0, 54, 0, 0, 0};
	UC bmpinfoheader[40] = {40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 24, 0};
	UC bmppad[3] = {0, 0, 0};
	bmpfileheader[ 2] = (UC)(filesize    );
	bmpfileheader[ 3] = (UC)(filesize>> 8);
	bmpfileheader[ 4] = (UC)(filesize>>16);
	bmpfileheader[ 5] = (UC)(filesize>>24);
	bmpinfoheader[ 4] = (UC)(       w    );
	bmpinfoheader[ 5] = (UC)(       w>> 8);
	bmpinfoheader[ 6] = (UC)(       w>>16);
	bmpinfoheader[ 7] = (UC)(       w>>24);
	bmpinfoheader[ 8] = (UC)(       h    );
	bmpinfoheader[ 9] = (UC)(       h>> 8);
	bmpinfoheader[10] = (UC)(       h>>16);
	bmpinfoheader[11] = (UC)(       h>>24);
	//writing file
	FILE* f = fopen(outfile,"wb");
	fwrite(bmpfileheader, sizeof(UC), 14, f);
	fwrite(bmpinfoheader, sizeof(UC), 40, f);
	for(int i=0; i<h; i++){
		fwrite(img + (w * (h-i-1) * 3), 3 , w, f);
		fwrite(bmppad, 1, (4 - (w*3)%4) % 4, f);
	}
	free(img);
	fclose(f);
}

void main(int argc, char **argv){
	if(argc == 1){
		printf("No Input Filename.\n");
		return;
	}
	FILE* fp = fopen(argv[1], "rb");
	
	UC marker[2];
	while(fread(marker, sizeof(UC), 2, fp) == 2){
		if(marker[0] == (UC)0xFF){ 
			switch(marker[1]){
				case (UC)0xDB:{ //DQT
					UC buffer[128];
					fread(buffer, sizeof(UC), 2, fp); //Lq
					int Lq = buffer[0] * 256 + buffer[1];
					int num_table = (Lq-2) / 65;
					
					//read DQT
					for(int i=0; i<num_table; i++){
						fread(buffer, sizeof(UC), 1, fp); //Pq & Tq
						int Pq = (buffer[0] & 0xF0) >> 4; //precision
						int Tq = buffer[0] & 0x0F; //table id
						int precision = (Pq==0) ? 1 : 2; //1 or 2 byte
						
						fread(buffer, sizeof(UC), 64 * precision, fp);
						for(int pos=0; pos<64; pos++)
							JPEG.DQT[Tq][pos] = (precision==1) ? buffer[pos] : buffer[pos*2]*256 + buffer[pos*2]; 
						//do zigzag after dequantization
					}
					break;
				}
				case (UC)0xC0:{ //SOF0
					UC buffer[2];
					fread(buffer, sizeof(UC), 2, fp); //Lf
					fread(buffer, sizeof(UC), 1, fp); //P, assume baseline = 0
					fread(buffer, sizeof(UC), 2, fp); //Y
					JPEG.Y = buffer[0] * 256 + buffer[1];
					fread(buffer, sizeof(UC), 2, fp); //X
					JPEG.X = buffer[0] * 256 + buffer[1];
					fread(buffer, sizeof(UC), 1, fp); //Nf, assume 3
					JPEG.Nf = buffer[0];
					
					JPEG.Hmax = 0;
					JPEG.Vmax = 0;
					for(int i=0; i<JPEG.Nf; i++){
						fread(buffer, sizeof(UC), 1, fp); //Ci, channel id
						int Ci = buffer[0] - 1; //shift into (0~2)
						fread(buffer, sizeof(UC), 1, fp); //H & V
						JPEG.H[Ci] = (buffer[0] & 0xF0) >> 4;
						JPEG.V[Ci] = buffer[0] & 0x0F;
						if(JPEG.H[Ci] > JPEG.Hmax)
							JPEG.Hmax = JPEG.H[Ci];
						if(JPEG.V[Ci] > JPEG.Vmax)
							JPEG.Vmax = JPEG.V[Ci];
						fread(buffer, sizeof(UC), 1, fp); //Tqi
						JPEG.Tqi[Ci] = buffer[0];
					}
					break;
				}
				case (UC)0xC4:{ //DHT
					UC buffer[2];
					fread(buffer, sizeof(UC), 2, fp); //Lh
					int Lh = buffer[0] * 256 + buffer[1];
					Lh -= 2;
					while(Lh > 0){
						fread(buffer, sizeof(UC), 1, fp); //Tc & Th
						int type = (buffer[0] & 0xF0) >> 4; //Tc
						int id = buffer[0] & 0x0F; //Th
						int V[16]; //numbers of various length codewords
						int sum = 0; //the sum of V 
						for(int i=0; i<16; i++){ //read 16 codewords
							fread(buffer, sizeof(UC), 1, fp); //Vi
							V[i] = (int)buffer[0];
							sum += V[i];
						}	

						UC* ptr = malloc( sum * sizeof(UC) ); //malloc enough space to load Huffman table
						fread(ptr, sizeof(UC), sum, fp);
				
						JPEG.DHT[id][type] = new_node();
						
						int len = 0; //the len of previous V[i]
						int code = 0; //encoded code
						for(int i=0; i<16; i++){ 
							ptr += len;
							len = V[i];
							for(int j=0; j<len; j++){
								int value = (int)(*(ptr + j)); //current value
								insert_node(JPEG.DHT[id][type], code, i+1, value); //insert value into Huffman tree
								code += 1;
							}
							code = code << 1;
						}
						Lh -= (17 + sum); //minus the length with the length already read
					}
					break;
				}
				case (UC)0xDA:{ //SOS
					UC buffer[3];
					fread(buffer, sizeof(UC), 2, fp); //Ls
					fread(buffer, sizeof(UC), 1, fp); // Ns, should be 3
					int Ns = (int)buffer[0];
					
					for(int i=0; i<Ns; i++){
						fread(buffer, sizeof(UC), 1, fp); //Csi, channel id
						int Csi = buffer[0] - 1; //shift into (0~2)
						fread(buffer, sizeof(UC), 1, fp); //Tdi & Tai
						JPEG.Tdi[Csi] = (buffer[0] & 0xF0) >> 4;
						JPEG.Tai[Csi] = buffer[0] & 0x0F;
					}
					fread(buffer, sizeof(UC), 3, fp); //three bytes discarded

					scan(fp); //scan compressed data
					break;
				}
			}
		}
		fseek(fp, -1, SEEK_SET);
	}

	YCbCr_to_RGB(); //convert to RGB
	write_bmp(argv[1]); //write bmp
	
	return;
}
