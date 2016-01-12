
/*!
 ***********************************************************************
 *  \file
 *     decoder_test.c
 *  \brief
 *     H.264/AVC decoder test 
 *  \author
 *     Main contributors (see contributors.h for copyright, address and affiliation details)
 *     - Yuwen He       <yhe@dolby.com>
 ***********************************************************************
 */
#include <sys/stat.h>

#include "win32.h"
#include "h264decoder.h"
#include "configfile.h"

long NumberOfMV=0;
long NumberOfFrame=0;
long NumberOfBMV=0;
long NumberOfPMV=0;

char g_encrypt_file[FILE_NAME_SIZE] = "../vediofile/decoder/encrypt.ec";
char g_key_log_file[FILE_NAME_SIZE] = "../vediofile/decoder/key_log.ec";

FILE * g_encrypt_fileh;

#if H264_KEY_LOG
FILE * g_key_log_fileh;
#endif


static void Configure(InputParameters *p_Inp, int ac, char *av[])
{
  memset(p_Inp, 0, sizeof(InputParameters));
  
  ParseCommand(p_Inp, ac, av);

  fprintf(stdout,"----------------------------- JM %s %s -----------------------------\n", VERSION, EXT_VERSION);
  if(!p_Inp->bDisplayDecParams)
  {
    fprintf(stdout,"--------------------------------------------------------------------------\n");
    fprintf(stdout," Input H.264 bitstream                  : %s \n",p_Inp->infile);
    fprintf(stdout,"--------------------------------------------------------------------------\n");
  }
  
}

/*!
 ***********************************************************************
 * \brief
 *    main function for JM decoder
 ***********************************************************************
 */
int main(int argc, char **argv)
{
	struct timeval start, end;
	gettimeofday( &start, NULL );
	
  int iRet;
  //DecodedPicList *pDecPicList;
  int iFramesDecoded=0;
  InputParameters InputParams;

  if((g_encrypt_fileh = fopen(g_encrypt_file, "w+")) == NULL)
  {
	  printf("encrypt file open fail!\n");
	  return -1;
  }

#if H264_KEY_LOG 
  if((g_key_log_fileh = fopen(g_key_log_file, "a+")) == NULL)
  {
	  printf("encrypt file open fail!\n");
	  return -1;
  }
#endif  

  //init_time();

  //get input parameters;
  Configure(&InputParams, argc, argv);

  //open decoder;
  iRet = OpenDecoder(&InputParams);	//打开trace,h264bitstream文件
  if(iRet != DEC_OPEN_NOERR)
  {
    fprintf(stderr, "Open encoder failed: 0x%x!\n", iRet);
    return -1; //failed;
  }

	p_Dec->nalu_pos_array = calloc(200,sizeof(int));
	
  do
  {
    iRet = DecodeOneFrame();	
  }//while((iRet == DEC_SUCCEED) && ((p_Dec->p_Inp->iDecFrmNum==0) || (iFramesDecoded<p_Dec->p_Inp->iDecFrmNum)));
  while(iRet == DEC_SUCCEED);
	
  iRet = FinitDecoder();
  iRet = CloseDecoder();	//包含report输出

  
  printf("%d frames are decoded.\n", iFramesDecoded);
  printf("%ld MVs found!\n", NumberOfMV);
  printf("%ld P MVs found!\n", NumberOfBMV);
  printf("%ld B MVs found!\n", NumberOfPMV);
  
#if TRACE
	printf("defined trace!\n");
#endif  

#if H264_KEY_LOG
	char s[200];
	snprintf(s,200,"max_mvd_BOffset: %4d, min_KeyData_Len: %2d, max_KeyData_Len: %4d\n",
	p_Dec->max_MVD_BOffset,p_Dec->min_KeyData_Len,p_Dec->max_KeyData_Len);		
	//fwrite(s,strlen(s),1,g_key_log_fileh);

	char s2[200];

	int i =0;
	for(;i<p_Dec->nalu_pos_array_idx;++i)
	{
		snprintf(s2,200,"%2d: %5d\n",i,p_Dec->nalu_pos_array[i]);		
		fwrite(s2,strlen(s2),1,g_key_log_fileh);		
	}
	fclose(g_key_log_fileh);	
#else	
	fclose(g_encrypt_fileh);
#endif

	gettimeofday( &end, NULL );
	long int time_us = 1000000 * ( end.tv_sec - start.tv_sec ) + end.tv_usec - start.tv_usec;
	printf("run time: %ld us\n",time_us);
  return 0;
}


