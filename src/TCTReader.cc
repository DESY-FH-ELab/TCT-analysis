/**
 * \file
 * \brief Implementation of the TCTReader and TCTWaveform classes. © Particulars d.o.o.
 * \details The code below was originally written by © Particulars d.o.o. and was called TCTAnalyse.
 * If needed one can download the original version from particulars.si.
 *
 */


#include "TCTReader.h"
#include "TMath.h"
#include "TPaveText.h"

//ClassImp(TCTReader);

TCTReader::TCTReader(char *FileNameInp, Float_t time0, Int_t Bin)
{
    FileName = FileNameInp;
    ////////////////////////////////////////////////////////////////////////////
    // Class for manipulation of the 3D/Position Sensitive TCT Measurements
    //
    // char *FileName;  name of the file with data
    // Float_t time0 ;  time shift of the points. Used mainly to set the arrival of the laser pulse at t=0 ns;
    // Int Bin; Selects the timeformat of the measurements;
    //          0 - ascii (default) older format has a type 11 while newer has type 22
    //          1 - binary (prefered in new measurements)
    //          2 - binary (little endian)
    // Example of use :
    // // Convert into TCTWaveform along projection
    // TCTReader aa("../Meritve/scanz-grobo-1.tct", 92.2,1); // The second parameter is to set the scale such that signal start at t=0;
    // aa.CorrectBaseLine();   // Baseline correction
    // aa.PrintInfo();         // Information about the read data
    // // example of projection the data along Y, 41 point at indexes of Z=23,X=0, U1=0, U2=0;
    // TCTWaveform *wf=aa.Projection(0,1,0,0,23,0,0,41);
    // strcpy(wf.suffix," #mum");
    // strcpy(wf.prefix,"Y=");
    // //Draw waveforms and substract the 0 waveform from the rest (cancel oscilations
    // wf.DrawMulti(-1,40,1,40,4,0,0);

    histo1 = NULL;
    histo2 = NULL;
    histo3 = NULL;
    histo4 = NULL;

    Int_t i,j,Cs,Us,Ss,ofs=0;
    Char_t filef[5];
    float header[100];
    for(i=0;i<4;i++) WFOnOff[i]=0;
    Date=TArrayI(6);
    User=NULL;
    Comment=NULL;
    Sample=NULL;

    if(Bin==2) BLE_CODE=0; else BLE_CODE=1;
    // BLE_CODE=0;
    if(!Bin) sprintf(filef,"r+"); else sprintf(filef,"rb+");
    if((in=fopen((const Char_t *)FileName,filef))==NULL) {printf("\n Error opening file for reading\n"); return;}

    if(!Bin)  // read ASCII file
    {
        fscanf(in,"%d",&type); // check the file type
        if(!(type==11 || type==22 || type==33 ))  // if it is something else exit
        {
            printf("Can not read other formats than waveform: %d!\n",type);
        }
        else
        {
            fscanf(in,"%d %d %d %d %d %d\n",&Date[0],&Date[1],&Date[2],&Date[3],&Date[4],&Date[5]);
            fscanf(in,"%d\n",&abstime);
            fscanf(in,"%f %f %d\n",&x0,&dx,&Nx);
            fscanf(in,"%f %f %d\n",&y0,&dy,&Ny);
            fscanf(in,"%f %f %d\n",&z0,&dz,&Nz);

            if(type==22) fscanf(in,"%d %d %d\n",&WFOnOff[0],&WFOnOff[1],&WFOnOff[2]); //Read in the wafeform off on
            if(type==33) fscanf(in,"%d %d %d %d\n",&WFOnOff[0],&WFOnOff[1],&WFOnOff[2],&WFOnOff[3]); //Read in the wafeform off on

            fscanf(in,"%d",&NU1); U1=TArrayF(NU1);
            for(i=0;i<NU1;i++) fscanf(in,"%f",&U1[i]); // printf("%d %d %f\n",NU1,i,U1[i]);}

            if(type!=11) fscanf(in,"%d",&NU2); else NU2=1;
            U2=TArrayF(NU2);
            if(type!=11) for(i=0;i<NU2;i++) fscanf(in,"%f",&U2[i]);

            I2=TArrayF(NU2*NU1);
            I1=TArrayF(NU2*NU1);

            fscanf(in,"%f %f %d\n",&t0,&dt,&NP);

            numxyz=Nx*Ny*Nz;

            for(i=0;i<8;i++) xyz[i]=new Float_t [numxyz*NU1*NU2];

            if(WFOnOff[0]) {histo1 =new TClonesArray("TH1F",numxyz*NU1*NU2); histo1->BypassStreamer(kFALSE);}
            if(WFOnOff[1]) {histo2 =new TClonesArray("TH1F",numxyz*NU1*NU2); histo2->BypassStreamer(kFALSE);}
            if(WFOnOff[2]) {histo3 =new TClonesArray("TH1F",numxyz*NU1*NU2); histo3->BypassStreamer(kFALSE);}
            if(WFOnOff[3]) {histo4 =new TClonesArray("TH1F",numxyz*NU1*NU2); histo4->BypassStreamer(kFALSE);}
        }
        ReadWFs(time0);
    }
    else
    {

        int read=fread((void *)header,sizeof(float),1,in);
        if(BLE_CODE) swooip(header,read);

        //get the file type
        type=(Int_t) header[0];
        std::cout<<"File type = "<<type<<" - reading \n"; // print the file type

        //read in the buffer
        rewind(in);
        read=fread((void *)header,sizeof(float),100,in);

        //     for(int ii=0;ii<50;ii++) printf("%d %f\n",ii,header[ii]);
        if(BLE_CODE) swooip(header,read);
        //get date
        for(i=0;i<6;i++) Date[i]=(Int_t) header[i+1];

        //get absolute time
        abstime=(int) header[7];
        //get moving matrix
        x0=header[8];  dx=header[9]; Nx=(int)header[10];
        y0=header[11]; dy=header[12]; Ny=(int)header[13];
        z0=header[14]; dz=header[15]; Nz=(int)header[16];

        //adjust the reading of the header
        if(type==33) ofs=1; else ofs=0;
        //get the on/off channels
        for(i=0;i<3+ofs;i++) WFOnOff[i]=(int)header[17+i];
        //get number of voltage steps for the first power source
        //     printf("%f %d\n",header[20+ofs],header[20+ofs]);
        NU1=(int)header[20+ofs]; U1=TArrayF(NU1);
        //get voltage steps for first power source
        for(i=0;i<NU1;i++) U1[i]=header[21+ofs+i];
        //get number of voltage steps for the second power source
        NU2=(int)header[21+NU1+ofs];U2=TArrayF(NU2);
        //get voltage steps for first power source
        for(i=0;i<NU2;i++) U2[i]=header[22+ofs+NU1+i];
        //get time scale
        t0=header[22+ofs+NU1+NU2]; if(TMath::Abs(t0)>1e-3) t0*=1e-9;
        dt=header[23+ofs+NU1+NU2]; if(TMath::Abs(dt)>1e-3) dt*=1e-9;
        NP=(int)header[24+ofs+NU1+NU2];
        //////    Header information coded from type=30 on
        switch(type)
        {
        case 33:
            T=header[25+ofs+NU1+NU2];
            Source=(int)header[26+ofs+NU1+NU2];
            //rewind to the appropriate position
            fseek(in,(27+NU1+NU2+ofs)*sizeof(Float_t),SEEK_SET);
            fread(&Us,sizeof(int),1,in); if(BLE_CODE) swooip((float *) &Us,1);
            User=new Char_t[Us+1];
            fread(User,sizeof(Char_t),Us,in);
            User[Us]='\0';
            fread(&Ss,sizeof(int),1,in); if(BLE_CODE)  swooip((float *) &Ss,1);
            Sample=new Char_t[Ss+1];
            fread(Sample,sizeof(Char_t),Ss,in);
            Sample[Ss]='\0';
            fread(&Cs,sizeof(int),1,in); if(BLE_CODE)  swooip((float *) &Cs,1);
            Comment=new Char_t[Cs+1];
            fread(Comment,1,Cs,in);
            Comment[Cs]='\0';
            //     fseek(in,(28+ofs+NU1+NU2+Us+Ss+Cs)*sizeof(Float_t),SEEK_SET);
            break;
        case 22:
            fseek(in,(25+NU1+NU2+ofs)*sizeof(Float_t),SEEK_SET);
            break;
        }

        ////////////////////////////////////////////
        //initializa current arrays
        I2=TArrayF(NU2*NU1);
        I1=TArrayF(NU2*NU1);
        //number of steps
        numxyz=Nx*Ny*Nz;
        for(i=0;i<8;i++) xyz[i]=new Float_t [numxyz*NU1*NU2];
        //intitialize histograms
        if(WFOnOff[0]) {histo1 =new TClonesArray("TH1F",numxyz*NU1*NU2); histo1->BypassStreamer(kFALSE);}
        if(WFOnOff[1]) {histo2 =new TClonesArray("TH1F",numxyz*NU1*NU2); histo2->BypassStreamer(kFALSE);}
        if(WFOnOff[2]) {histo3 =new TClonesArray("TH1F",numxyz*NU1*NU2); histo3->BypassStreamer(kFALSE);}
        if(WFOnOff[3]) {histo4 =new TClonesArray("TH1F",numxyz*NU1*NU2); histo4->BypassStreamer(kFALSE);}
        //for(i=0;i<50;i++) printf("%d %f\n",i,header[i]);
        ReadWFsBin(time0);

    }

    RefInd=-1;
    //Setting the color map
}

TCTReader::~TCTReader()
{
    if(histo1) {
        delete histo1;
        histo1 = NULL;
    }
    if(histo2) {
        delete histo2;
        histo2 = NULL;
    }
    if(histo3) {
        delete histo3;
        histo3 = NULL;
    }
    if(histo4) {
        delete histo4;
        histo4 = NULL;
    }
    delete User;
    delete Sample;
    delete Comment;
    for(int i=0;i<8;i++) delete xyz[i];
}

void  TCTReader::ReadWFsBin(Float_t time0)
{
    // read in binary waveforms
    Int_t i,ii,j,k,q,r,numread;
    Float_t data,tU1,tU2,tI1,tI2;
    Char_t hisname1[100];
    Char_t hisname2[100];
    Char_t hisname3[100];
    Char_t hisname4[100];
    TClonesArray &entryp1 = *histo1;
    TClonesArray &entryp2 = *histo2;
    TClonesArray &entryp3 = *histo3;
    TClonesArray &entryp4 = *histo4;

    Float_t buf[10000];

    for(q=0;q<NU1;q++)
    {
        for(r=0;r<NU2;r++)
        {

            fread((void *)buf,sizeof(Float_t),4,in); if(BLE_CODE) swooip(buf,4);
            tU1=buf[0]; tU2=buf[1]; tI1=buf[2]; tI2=buf[3];
            U1[q]=tU1; I1[r+q*NU2]=tI1;
            U2[r]=tU2; I2[r+q*NU2]=tI2;

            //  printf("%d %d :: %f %f %f %f\n",r,q,tU1,tU2,tI1,tI2);
            for(i=0;i<numxyz;i++)
            {

                ii=i+numxyz*r+(NU2*numxyz)*q;

                fread((void *)buf,sizeof(Float_t),4,in);  if(BLE_CODE) swooip(buf,4);
                // printf("%d :: %f %f %f %f\n",ii,buf[0],buf[1],buf[2],buf[3]);

                for(j=0;j<4;j++) xyz[j][ii]=buf[j];
                xyz[7][ii]=xyz[3][ii];
                xyz[3][ii]=tU1; xyz[4][ii]=tU2;
                xyz[5][ii]=tI1; xyz[6][ii]=tI2;


                if(WFOnOff[0])
                {
                    sprintf(hisname1,"Ch. 1:x=%.6e,y=%.6e,z=%.6e,U1=%4.2f, U2=%4.2f ",xyz[0][ii],xyz[1][ii],xyz[2][ii],xyz[3][ii],xyz[4][ii]);
                    new(entryp1[ii]) TH1F((const Char_t *)(hisname1),(const Char_t *)(hisname1),NP,t0*1e9-time0,(NP*dt+t0)*1e9-time0);
                }

                if(WFOnOff[1])
                {
                    sprintf(hisname2,"Ch. 2:x=%.6e,y=%.6e,z=%.6e,U1=%4.2f, U2=%4.2f ",xyz[0][ii],xyz[1][ii],xyz[2][ii],xyz[3][ii],xyz[4][ii]);
                    new(entryp2[ii]) TH1F((const Char_t *)(hisname2),(const Char_t *)(hisname2),NP,t0*1e9-time0,(NP*dt+t0)*1e9-time0);
                }

                if(WFOnOff[2])
                {
                    sprintf(hisname3,"Ch. 3:x=%.6e,y=%.6e,z=%.6e,U1=%4.2f, U2=%4.2f ",xyz[0][ii],xyz[1][ii],xyz[2][ii],xyz[3][ii],xyz[4][ii]);
                    new(entryp3[ii]) TH1F((const Char_t *)(hisname3),(const Char_t *)(hisname3),NP,t0*1e9-time0,(NP*dt+t0)*1e9-time0);
                }

                if(WFOnOff[3])
                {
                    sprintf(hisname4,"Ch. 4:x=%.6e,y=%.6e,z=%.6e,U1=%4.2f, U2=%4.2f ",xyz[0][ii],xyz[1][ii],xyz[2][ii],xyz[3][ii],xyz[4][ii]);
                    new(entryp4[ii]) TH1F((const Char_t *)(hisname4),(const Char_t *)(hisname4),NP,t0*1e9-time0,(NP*dt+t0)*1e9-time0);
                }


                for(k=0;k<4;k++)
                {
                    if(WFOnOff[k]==1)
                    {
                        //	printf("reading ...... %d %d ... ",k,ftell(in));
                        numread=fread(buf,sizeof(Float_t),NP,in);  if(BLE_CODE) swooip(buf, NP); // printf("%f %f %f\n",buf[NP-3],buf[NP-2],buf[NP-1]);
                        //	printf("read ...... %d(%d)\n",numread,ii);  if(ii==17) return;
                        //fread(buf,sizeof(Float_t),NP,in);  swooip(buf, NP);
                        for(j=0;j<NP;j++)
                        {
                            //if(ii<2) if(j%20!=0) printf("%5.2f", buf[j]); else printf("%5.2f\n", buf[j]);
                            switch(k)
                            {
                            case 0:
                                ((TH1F*)entryp1[ii])->SetBinContent(j+1,buf[j]);
                                //if(j==0) printf("ii=%d k=%d j=%d, buf=%f\n",ii,k,j,buf[j]*1e-3);
                                //if(ii=0) printf("ii=%d k=%d j=%d, buf=%f\n",ii,k,j,buf[j]*1e-3);
                                break;
                            case 1:
                                ((TH1F*)entryp2[ii])->SetBinContent(j+1,buf[j]);
                                break;
                            case 2:
                                ((TH1F*)entryp3[ii])->SetBinContent(j+1,buf[j]);
                                break;
                            case 3:
                                ((TH1F*)entryp4[ii])->SetBinContent(j+1,buf[j]);
                                break;
                            }
                        }

                        //		printf("k=%d \n",k);

                    }

                }

            }

        }

    }

    //fclose(in);
}



void  TCTReader::ReadWFs(Float_t time0)
{
    Int_t i,ii,j,k,q,r;
    Float_t data,tU1,tU2,tI1,tI2;
    Char_t hisname1[100];
    Char_t hisname2[100];
    Char_t hisname3[100];
    Char_t hisname4[100];
    TClonesArray &entryp1 = *histo1;
    TClonesArray &entryp2 = *histo2;
    TClonesArray &entryp3 = *histo3;
    TClonesArray &entryp4 = *histo4;


    for(q=0;q<NU1;q++)
    {
        for(r=0;r<NU2;r++)
        {

            if(type!=11)
            {
                fscanf(in,"%f %f %f %f",&tU1,&tU2,&tI1,&tI2);
                // printf("%d %d :: %f %f %f %f\n",r,q,tU1,tU2,tI1,tI2);
                U1[q]=tU1; I1[r+q*NU2]=tI1;
                U2[r]=tU2; I2[r+q*NU2]=tI2;
            }


            for(i=0;i<numxyz;i++)
            {

                ii=i+numxyz*r+(NU2*numxyz)*q;

                if(type==22 || type==33)
                {
                    for(j=0;j<4;j++) fscanf(in,"%f",&xyz[j][ii]); xyz[7][ii]=xyz[3][ii];
                    xyz[3][ii]=tU1; xyz[4][ii]=tU2; xyz[5][ii]=tI1; xyz[6][ii]=tI2;
                }
                if(type==11) for(j=0;j<5;j++) fscanf(in,"%f",&xyz[j][ii]);

                if(WFOnOff[0])
                {
                    sprintf(hisname1,"Ch. 1:x=%.6e,y=%.6e,z=%.6e,U1=%4.2f, U2=%4.2f ",xyz[0][ii],xyz[1][ii],xyz[2][ii],xyz[3][ii],xyz[4][ii]);
                    new(entryp1[ii]) TH1F((const Char_t *)(hisname1),(const Char_t *)(hisname1),NP,t0*1e9-time0,(NP*dt+t0)*1e9-time0);
                }
                if(WFOnOff[1])
                {
                    sprintf(hisname2,"Ch. 2:x=%.6e,y=%.6e,z=%.6e,U1=%4.2f, U2=%4.2f ",xyz[0][ii],xyz[1][ii],xyz[2][ii],xyz[3][ii],xyz[4][ii]);
                    new(entryp2[ii]) TH1F((const Char_t *)(hisname2),(const Char_t *)(hisname2),NP,t0*1e9-time0,(NP*dt+t0)*1e9-time0);
                }
                if(WFOnOff[2])
                {
                    sprintf(hisname3,"Ch. 3:x=%.6e,y=%.6e,z=%.6e,U1=%4.2f, U2=%4.2f ",xyz[0][ii],xyz[1][ii],xyz[2][ii],xyz[3][ii],xyz[4][ii]);
                    new(entryp3[ii]) TH1F((const Char_t *)(hisname3),(const Char_t *)(hisname3),NP,t0*1e9-time0,(NP*dt+t0)*1e9-time0);
                }

                if(WFOnOff[3])
                {
                    sprintf(hisname4,"Ch. 4:x=%.6e,y=%.6e,z=%.6e,U1=%4.2f, U2=%4.2f ",xyz[0][ii],xyz[1][ii],xyz[2][ii],xyz[3][ii],xyz[4][ii]);
                    new(entryp4[ii]) TH1F((const Char_t *)(hisname4),(const Char_t *)(hisname4),NP,t0*1e9-time0,(NP*dt+t0)*1e9-time0);
                }


                for(k=0;k<4;k++)
                {
                    if(WFOnOff[k]==1)
                    {
                        for(j=0;j<NP;j++)
                        {
                            fscanf(in,"%e",&data); //printf("%e ",data);
                            switch(k)
                            {
                            case 0:
                                ((TH1F*)entryp1[ii])->SetBinContent(j+1,data);
                                break;
                            case 1:
                                ((TH1F*)entryp2[ii])->SetBinContent(j+1,data);
                                break;
                            case 2:
                                ((TH1F*)entryp3[ii])->SetBinContent(j+1,data);
                                break;
                            case 3:
                                ((TH1F*)entryp4[ii])->SetBinContent(j+1,data);
                                break;

                            }
                        }
                        //	  printf("k=%d \n",k);
                    }
                }
            }
        }
    }
    //fclose(in);
}

TH1F *TCTReader::GetHA(int ch , int index)
{
    TH1F *his;
    switch(ch)
    {
    case 0:
        his=(TH1F *)histo1->At(index);
        his->SetLineColor(1);
        break;
    case 1: his=(TH1F *)histo2->At(index);
        his->SetLineColor(2);
        break;
    case 2: his=(TH1F *)histo3->At(index);
        his->SetLineColor(4);
        break;
    case 3: his=(TH1F *)histo4->At(index);
        his->SetLineColor(5);
        break;
    default: his=NULL; break;
    }
    his->GetXaxis()->SetTitle("t [ns]");
    his->GetYaxis()->SetTitle("I [arb.]");

    return(his);
}

TH1F *TCTReader::GetHA(Int_t ch , Int_t x, Int_t y, Int_t z, Int_t nu1, Int_t nu2)
{
    return(GetHA(ch,indx(x,y,z,nu1,nu2)));
}


Int_t TCTReader::indx(int x, int y, int z, int nu1, int nu2)
{
    // Function returns the index of the position of the waveform corresponding to x,y,z,nu1,nu2 in the linear array of waveforms.
    // int x; index in the X direction
    // int y; index in the Y direction
    // int z; index in the Z direction
    // int nu1; index of desired votlage (1)
    // int nu2; index of desired voltage (2)
    if(x>Nx-1 || x<0) {printf("index x: out of range\n"); return 0;}
    if(y>Ny-1 || y<0) {printf("index y: out of range\n"); return 0;}
    if(z>Nz-1 || z<0) {printf("index z: out of range\n"); return 0;}
    if(nu1>NU1-1 || nu1<0) {printf("index nu1: out of range\n"); return 0;}
    if(nu2>NU2-1 || nu2<0) {printf("index nu2: out of range\n"); return 0;}

    return( (x+Nx*y+(Nx*Ny)*z)+numxyz*nu2+(NU2*numxyz)*nu1 );
};

void TCTReader::cords(int x, int y, int z, int nu1, int nu2)
{
    Int_t ix=indx(x,y,z,nu1,nu2);
    printf("X(%d)=%f, Y(%d)=%f, Z(%d)=%f :: ",ix,xyz[0][ix],ix,xyz[1][ix],ix,xyz[2][ix]);
    printf("%d-> U1=%f, U2=%f, I1=%f, I2=%f :: %4.0f\n",ix,xyz[3][ix],xyz[4][ix],xyz[5][ix],xyz[6][ix],xyz[7][ix]);
};

Float_t TCTReader::GetWidth(TH1F *his, Int_t &lbin, Int_t &hbin, Float_t tleft, Float_t tright, Float_t minwidth, Float_t maxwidth, Int_t percentage)
{
    //his -> histogram
    //tleft -> left margin, default 0
    //tright -> right margin, default 25
    //minwidth -> min signal width, default -1111
    //maxwidth -> max signal width, default -1111
    //percentage -> height of min and max, default 50 (50%)

    //return Dt -> ok
    //return -No of points -> if minwidth & maxwidth are default values
    //return -53 -> error

    Float_t q [100];
    Int_t j = 0,i,i1,i2,index1,index2,pp;
    Float_t val,valn;
    Float_t Dt=0,pDt;
    Float_t max=0,pmax,max2;
    //Float_t *res=new Float_t [3];

    for(i=his->GetXaxis()->FindBin(tleft);i<his->GetXaxis()->FindBin(tright);i++)
    {
        pmax=his->GetBinContent(i);

        if(max<pmax)
        {
            max=pmax;
            max2=max*percentage/100;
        }
    }


    for(i=his->GetXaxis()->FindBin(tleft);i<his->GetXaxis()->FindBin(tright);i++)

    {
        valn=his->GetBinContent(i+1);
        val=his->GetBinContent(i);

        if((val<=max2 && valn>max2) || (val>max2 && valn<=max2))
        {
            q[j]=his->GetBinCenter(i);
            j++;
        }
    }

    //	printf("Maximum :: %f\n",max);

    if(j==2)
    {
        //printf("Interval %f %f\n",q[0],q[1]);
        lbin=his->GetXaxis()->FindBin(q[0]);
        hbin=his->GetXaxis()->FindBin(q[1]);
        return Dt = q[1]-q[0];

    }

    else if(j>2)
    {
        if((minwidth==-1111 || maxwidth==-1111)||(minwidth==0 || maxwidth==0)) return -j;
        else
        {
            index1=j-1;
            while(index1>0)
            {
                index2 = index1-1;
                while(index2>=0)
                {
                    pDt=q[index1]-q[index2];
                    //printf("q[%d] :: %f, q[%d] :: %f, pDt :: %f\n",index1,q[index1],index2,q[index2],pDt);
                    if(maxwidth>=pDt && pDt>=minwidth)
                    {
                        if(pDt>Dt)
                        {
                            Dt=pDt;
                        }
                    }
                    index2--;
                }
                index1--;
            }
            //printf("Interval %d %d\n",index1,index2);
            lbin=his->GetXaxis()->FindBin(q[index1]);
            hbin=his->GetXaxis()->FindBin(q[index2]);
            return Dt;
        }
    }else return -53;
}


TH2F *TCTReader::Draw(Int_t ch, Int_t stype, Int_t mode, Int_t v1, Int_t v2, Int_t v3, Float_t tlow, Float_t thi)
{
    //  Int_t ch; channel number
    //  Int_t stype
    // case: 0 - select x,y
    // case: 1 - select x,z
    // case: 2 - select x,U1
    // case: 3 - select x,U2
    // case: 4 - select y,z
    // case: 5 - select y,U1
    // case: 6 - select y,U2
    // case: 7 - select z,U1
    // case: 8 - select z,U2
    // case: 9 - select U1,U2
    // Int_t mode;
    //          mode =0 -> Integral
    //          mode =1 -> maximum signal
    //          mode =2 -> minimum signal
    //          mode =3 -> width of the signal
    // Int_t v1,v2,v3; indexes of the remaining three parameters (order always follows x,y,z,U1,U2)
    // Float_t tlow, thi; Time windos of the plots
    //
    //
    // Example :
    // //plot channel 0 of the, yz plot, with indexesx=0, U1=0, U2=0, in the time windows [0,30 ns]
    //  TH2F *plot=aa.Draw(0,4,0,0,0,0,0,20);
    //  plot->Draw("COLZ");

    Char_t txt[100];
    Float_t integral,max,min,width;
    Int_t i,j,left,right,t,pX,pY,u1s,u1e,u2s,u2e;
    Int_t *s[5];
    Int_t lbin,hbin;

    TH2F *his2;
    TH1F *his;

    switch(stype)
    {
    case 0://x,y
        t = indx(0,0,v1,v2,v3);
        sprintf(txt,"2d his (Z[%d]=%f U1[%d]=%f U2[%d]=%f)",v1,xyz[2][t],v2,xyz[3][t],v3,xyz[4][t]);
        his2 = new TH2F(txt,txt,Nx,x0,Nx*dx+x0,Ny,y0,Ny*dy+y0);
        s[0] = &i,s[1] = &j,s[2] = &v1,s[3] = &v2,s[4] = &v3;
        pX=Nx; pY=Ny;
        break;
    case 1://x,z
        t = indx(0,v1,0,v2,v3);
        sprintf(txt,"2d his (Y[%d]=%f U1[%d]=%f U2[%d]=%f)",v1,xyz[1][t],v2,xyz[3][t],v3,xyz[4][t]);
        his2 = new TH2F(txt,txt,Nx,x0,Nx*dx+x0,Nz,z0,Nz*dz+z0);
        s[0] = &i,s[1] = &v1,s[2] = &j,s[3] = &v2,s[4] = &v3;
        pX=Nx; pY=Nz;
        break;
    case 2://x,U1
        t = indx(0,v1,v2,0,v3);
        sprintf(txt,"2d his (Y[%d]=%f Z[%d]=%f U2[%d]=%f)",v1,xyz[1][t],v2,xyz[2][t],v3,xyz[4][t]);
        if(U1[NU1-1]<U1[0]) {u1s=U1[NU1-1]; u1e=U1[0];} else { u1s=U1[0]; u1e=U1[NU1-1];}
        his2 = new TH2F(txt,txt,(Int_t) Nx,(Float_t) x0,Nx*dx+x0,NU1,u1s,u1e);
        s[0] = &i,s[1] = &v1,s[2] = &v2,s[3] = &j,s[4] = &v3;
        pX=Nx; pY=NU1;
        break;
    case 3://x,U2
        t = indx(0,v1,v2,v3,0);
        sprintf(txt,"2d his (Y[%d]=%f Z[%d]=%f U1[%d]=%f)",v1,xyz[1][t],v2,xyz[2][t],v3,xyz[3][t]);
        if(U2[NU2-1]<U2[0]) {u2s=U2[NU2-1]; u2e=U2[0];} else { u2s=U2[0]; u2e=U2[NU2-1];}
        his2 = new TH2F(txt,txt,Nx,x0,Nx*dx+x0,NU2,u2s,u2e);
        s[0] = &i,s[1] = &v1,s[2] = &v2,s[3] = &v3,s[4] = &j;
        pX=Nx; pY=NU2;
        break;
    case 4://y,z
        t = indx(v1,0,0,v2,v3);
        sprintf(txt,"2d his (X[%d]=%f U1[%d]=%f U2[%d]=%f)",v1,xyz[2][t],v2,xyz[3][t],v3,xyz[4][t]);
        his2 = new TH2F(txt,txt,Ny,y0,Ny*dy+y0,Nz,z0,Nz*dz+z0);
        s[0] = &v1,s[1] = &i,s[2] = &j,s[3] = &v2,s[4] = &v3;
        pX=Ny; pY=Nz;
        break;
    case 5://y,U1
        t = indx(v1,0,v2,0,v3);
        sprintf(txt,"2d his (X[%d]=%f Z[%d]=%f U2[%d]=%f)",v1,xyz[0][t],v2,xyz[2][t],v3,xyz[4][t]);
        if(U1[NU1-1]<U1[0]) {u1s=U1[NU1-1]; u1e=U1[0];} else { u1s=U1[0]; u1e=U1[NU1-1];}
        his2 = new TH2F(txt,txt,Ny,y0,Ny*dy+y0,NU1,u1s,u1e);
        s[0] = &v1,s[1] = &i,s[2] = &v2,s[3] = &j,s[4] = &v3;
        pX=Ny; pY=NU1;

        break;
    case 6://y,U2
        t = indx(v1,0,v2,v3,0);
        sprintf(txt,"2d his (X[%d]=%f Z[%d]=%f U1[%d]=%f)",v1,xyz[0][t],v2,xyz[2][t],v3,xyz[3][t]);
        if(U2[NU2-1]<U2[0]) {u2s=U2[NU2-1]; u2e=U2[0];} else { u2s=U2[0]; u2e=U2[NU2-1];}
        his2 = new TH2F(txt,txt,Ny,y0,Ny*dy+y0,NU2,u2s,u2e);
        s[0] = &v1,s[1] = &i,s[2] = &v2,s[3] = &v3,s[4] = &j;
        pX=Ny; pY=NU2;
        break;
    case 7://z,U1
        t = indx(v1,v2,0,0,v3);
        sprintf(txt,"2d his (X[%d]=%f Y[%d]=%f U2[%d]=%f)",v1,xyz[0][t],v2,xyz[1][t],v3,xyz[4][t]);
        if(U1[NU1-1]<U1[0]) {u1s=U1[NU1-1]; u1e=U1[0];} else { u1s=U1[0]; u1e=U1[NU1-1];}
        his2 = new TH2F(txt,txt,Nz,z0,Nz*dz+z0,NU1,u1s,u1e);
        s[0] = &v1,s[1] = &v2,s[2] = &i,s[3] = &j,s[4] = &v3;
        pX=Nz; pY=NU1;
        break;
    case 8://z,U2
        t = indx(v1,v2,0,v3,0);
        sprintf(txt,"2d his (X[%d]=%f Y[%d]=%f U1[%d]=%f)",v1,xyz[0][t],v2,xyz[1][t],v3,xyz[3][t]);
        if(U2[NU2-1]<U2[0]) {u2s=U2[NU2-1]; u2e=U2[0];} else { u2s=U2[0]; u2e=U2[NU2-1];}
        his2 = new TH2F(txt,txt,Nz,z0,Nz*dz+z0,NU2,u2s,u2e);
        s[0] = &v1,s[1] = &v2,s[2] = &i,s[3] = &v3,s[4] = &j;
        pX=Nz; pY=NU2;
        break;
    case 9://U1,U2
        t = indx(v1,v2,v3,0,0);
        sprintf(txt,"2d his (X[%d]=%f Y[%d]=%f Z[%d]=%f)",v1,xyz[0][t],v2,xyz[1][t],v3,xyz[2][t]);
        if(U2[NU2-1]<U2[0]) {u2s=U2[NU2-1]; u2e=U2[0];} else { u2s=U2[0]; u2e=U2[NU2-1];}
        if(U1[NU1-1]<U1[0]) {u1s=U1[NU1-1]; u1e=U1[0];} else { u1s=U1[0]; u1e=U1[NU1-1];}
        his2 = new TH2F(txt,txt,NU1,u1s,u1e,NU2,u2s,u2e);
        s[0] = &v1,s[1] = &v2,s[2] = &v3,s[3] = &i,s[4] = &j;
        pX=NU1; pY=NU2;
        break;
    }

    for(j=0;j<pY;j++)
        for(i=0;i<pX;i++)
        {

            switch(mode)
            {
            case 0:		// integral
                his=GetHA(ch,*s[0],*s[1],*s[2],*s[3],*s[4]);
                right=his->GetXaxis()->FindBin(thi);
                left=his->GetXaxis()->FindBin(tlow);
                integral=his->Integral(left,right);
                //printf("Vrednost integrala :: %f\n",integral);
                break;
            case 1:		// max

                his=GetHA(ch,*s[0],*s[1],*s[2],*s[3],*s[4]);
                integral=GetHA(ch,*s[0],*s[1],*s[2],*s[3],*s[4])->GetMaximum();
                width = GetWidth(his,lbin,hbin,0,40);
                //if(width > 0)
                //{
                //	printf("[i,j,width] = [%d  %d  %f]\n",i,j,width);
                //}
                break;
            case 2:		// min
                integral=GetHA(ch,*s[0],*s[1],*s[2],*s[3],*s[4])->GetMinimum();
                break;
            case 3:		// Dt
                his=GetHA(ch,*s[0],*s[1],*s[2],*s[3],*s[4]);
                integral = GetWidth(his,lbin,hbin,30,60,3,5);
                break;
            }
            //			   printf("i=%d,j=%d int=%f:: %d %d %d %d %d\n",i,j,integral,*s[0],*s[1],*s[2],*s[3],*s[4]);
            //printf("i=%d,j=%d int=%f:: %d %d %d %d %d\n",i,j,integral,*s[0],*s[1],*s[2],*s[3],*s[4]);
            his2->SetBinContent(i+1,j+1,integral);
        }

    max=his2->GetMaximum();
    min=his2->GetMinimum();

    if(TMath::Abs(max)>TMath::Abs(min))
        his2->SetMinimum(-TMath::Abs(max));
    else
        his2->SetMaximum(TMath::Abs(min));
    his2->Draw("SURF2");
    return his2;
}

TCTWaveform *TCTReader::Projection(Int_t num, Int_t *List)
{

    Float_t Delta;
    Int_t i=0,ix;
    TCTWaveform *MWF=new TCTWaveform(num);
    TH1F *his;
    for(i=0;i<num;i++)
    {
        ix=indx(List[i*6+1],List[i*6+2],List[i*6+3],List[i*6+4],List[i*6+5]);
        his=GetHA(List[i*6],ix);
        MWF->AddHisto(i,i,his);
    }
    strcpy(MWF->suffix,"");
    MWF->DrawMode=false;
    strcpy(MWF->prefix,"wf=");
    return MWF;
}

TCTWaveform *TCTReader::Projection(int ch, int dir,int x,int y,int z, int nu1, int nu2, int num)
{
    // Projection parameters
    // int ch  -> channel number
    // int dir -> direction of the projection
    // int x   -> x0 of the projection
    // int y   -> y0 of the projection
    // int z   -> z0 of the projection
    // int nu1 -> voltage 1
    // int nu1 -> voltage 2
    // int num -> number of wfs

    Float_t Delta;
    Int_t i=0,ix;
    TCTWaveform *MWF=new TCTWaveform(num);
    TH1F *his;
    for(i=0;i<num;i++)
    {
        switch(dir)
        {
        case 0: ix=indx(i+x,y,z,nu1,nu2); Delta=i*dx; strcpy(MWF->prefix,"x=");   strcpy(MWF->suffix," #mum"); break;
        case 1: ix=indx(x,y+i,z,nu1,nu2); Delta=i*dy; strcpy(MWF->prefix,"y=");   strcpy(MWF->suffix," #mum"); break;
        case 2: ix=indx(x,y,z+i,nu1,nu2); Delta=i*dz; strcpy(MWF->prefix,"z=");   strcpy(MWF->suffix," #mum"); break;
        case 3: ix=indx(x,y,z,nu1+i,nu2); Delta=U1[nu1+i]; strcpy(MWF->prefix,"U1="); strcpy(MWF->suffix," V");    break;
        case 4: ix=indx(x,y,z,nu1,nu2+i); Delta=U2[nu2+i]; strcpy(MWF->prefix,"U2="); strcpy(MWF->suffix," V");    break;
        default: ix=indx(i+x,y,z,nu1,nu2); Delta=i*dx; break;
        }
        his=GetHA(ch,ix);
        MWF->AddHisto(i,Delta,his);
    }
    MWF->DrawMode=false;
    //    strcpy(MWF->suffix," #mum");
    return MWF;
}

void TCTReader::DrawList(Int_t num, Int_t *List)
{
    Float_t Delta;
    Int_t i=0,ix;
    TH1F *his;
    for(i=0;i<num;i++)
    {
        his=GetHA(List[i*6],List[i*6+1],List[i*6+2],List[i*6+3],List[i*6+4],List[i*6+5]);
        if(i==0) his->DrawCopy(); else his->DrawCopy("SAME");
    }

}

void TCTReader::DrawList(Int_t num, Int_t *ListC,Int_t *ListP)
{
    // Function draws the list
    Float_t Delta;
    Int_t i=0,ix;
    TH1F *his;
    for(i=0;i<num;i++)
    {
        his=GetHA(ListC[i],ListP[i]);
        if(i==0) his->Draw(); else his->Draw("SAME");
    }

}

void TCTReader::CorrectBaseLine(Float_t xc)
{
    // Function corrects the baseline (DC offset) of all wafeforms
    // Float_t xc ; time denoting the start of the pulse
    //              correction factor is calculated from all the bins before xc
    Int_t right[4],left[4];
    Int_t i,j,k;
    Double_t sum=0,corr[4];
    TH1F *his[4];
    Int_t Num=numxyz*NU1*NU2; //number of all waveforms


    for(j=0;j<Num;j++)
    {
        if(j==0)  std::cout<<"Baseline correction ("<<Num<<" waveforms) :: ";

        if(WFOnOff[0]==1) his[0]=((TH1F *)histo1->At(j));
        if(WFOnOff[1]==1) his[1]=((TH1F *)histo2->At(j));
        if(WFOnOff[2]==1) his[2]=((TH1F *)histo3->At(j));
        if(WFOnOff[3]==1) his[3]=((TH1F *)histo4->At(j));

        for(i=0;i<4;i++)
        {
            if(WFOnOff[i]==1)
            {
                right[i]=his[i]->GetXaxis()->FindBin(xc);
                left[i]=1;
                his[i]->Integral(left[i],right[i]);
                corr[i]=his[i]->Integral(left[i],right[i])/(right[i]-left[i]);
            }
        }

        if(j%100==0) std::cout<<".";
        //printf("%d :: Baseline correction = %e , Integral before trigger=%e , Nbins=%d!\n",j,corr,his->Integral(left,right),right-left);

        for(k=0;k<3;k++)
            if(WFOnOff[k]==1)
                for(i=1;i<his[k]->GetNbinsX();i++)
                    his[k]->SetBinContent(i,his[k]->GetBinContent(i)-corr[k]);

    }

    std::cout<<" finished\n";

}


void TCTReader::PrintInfo()
{
    // Function prints the information about the class and its members

    Int_t i,j;
    printf("Format of the file %d\n",type);
    printf("*************************************\n");
    if(User!=NULL) printf("User: %s \n",User);
    if(Sample!=NULL) printf("Sample: %s \n",Sample);
    if(Comment!=NULL) printf("Comment: %s \n",Comment);
    printf("*************************************\n");
    printf("Date and time of the meaurement: %d.%d.%d %d:%d:%d\n",Date[0],Date[1],Date[2],Date[3],Date[4],Date[5]);
    printf("Active osciloscope ch: Ch1=%d Ch2=%d Ch3=%d Ch4=%d \n",WFOnOff[0],WFOnOff[1],WFOnOff[2],WFOnOff[3]);
    printf("Number of points %d (X=%d, Y=%d, Z=%d)\n",Nx*Ny*Nz,Nx,Ny,Nz);
    printf("Positions: r0=(%f,%f,%f) dr=(%f,%f,%f) \n",x0,y0,z0,dx,dy,dz);
    printf("Time scale: points=%d, t0=%e, dt=%e\n",NP,t0,dt);

    printf("Temperature: %f\n",T);
    printf("Type of generation: %4.0f\n",Source);

    printf("Voltages: NU1=%d , NU2=%d::\n",NU1,NU2);
    for(i=0;i<NU1;i++)
        for(j=0;j<NU2;j++)
            printf("U1,U2(%f,%f)::I1,I2(%e,%e)\n",U1[i],U2[j],I1[j+NU2*i],I2[j+NU2*i]); printf("\n");



}


void  TCTReader::swoo(char *a, char *b) {
    // byte swaping (LABVIEW,HPUX g++)<->(LINUX g++, WINNT cl)
    char c = *a;
    *a = *b; *b = c;
}

void TCTReader::swooip(float *in, int s) {
    // byte swaping (LABVIEW,HPUX g++)<->(LINUX g++, WINNT cl)
    char *sr, b;
    while(s--) {
        sr=(char *)in;
        swoo(&sr[0], &sr[3]);
        swoo(&sr[1], &sr[2]);
        in++;
    }
}

TCTWaveform::TCTWaveform(Int_t num)
{
    // Default Constructor
    // 		Int_t num    ; Number of aqusitions
    Multiple=num;
    histo =new TClonesArray("TH1F",Multiple);
    histo->BypassStreamer(kFALSE);
    Voltages=TArrayF(Multiple);
    Temperature=TArrayF(Multiple);
    Current=TArrayF(Multiple);
    Frequencies=TArrayF(Multiple);
    Date=TArrayF(6);
    Frequency=0;

    Date.Reset();
    Frequencies.Reset();
    Voltages.Reset();
    Date.Reset();
    Current.Reset();
    Temperature.Reset();
    DrawMode=true;
    strcpy(suffix," V");
    strcpy(prefix,"U=");

    pt=NULL;
}

TCTWaveform::~TCTWaveform()
{
    // Destructor
    if(pt!=NULL) delete pt;
    //Clear();
    if(histo) {
        delete histo;
        histo = NULL;
    }
}

void TCTWaveform::AddHisto(Int_t index,Float_t voltage,TH1F *hisin)
{
    // Add histogram (waveform)
    // 		Int_t Index     ;  Index of the waveform
    // 		Float_t Voltage ;  Voltage at Index
    // 		TH1F*           ;  pointer to histogram of the waveform
    // This function is internally used for reading the measurement files
    if(index>Multiple || index<0) printf("Array out of bonds\n"); else
    {
        Voltages[index]=voltage;
        TClonesArray &entryp = *histo;
        new(entryp[index]) TH1F();
        hisin->Copy(*entryp[index]);
    }
}

void TCTWaveform::AddHisto(Float_t voltage,TH1F *hisin)
{
    // Add histogram (waveform)
    // 	        Float_t Voltage ;  Voltage at Index
    // 		TH1F*           ;  pointer to histogram of the waveform
    // This function is internally used for reading the measurement files.
    // It is assumed that the array of waveforms is already initialized
    // with voltages specified
    Int_t index,i;
    for(i=0;i<Multiple;i++) if(Voltages[i]==voltage)
    {
        index=i;
        TClonesArray &entryp = *histo;
        new(entryp[index]) TH1F();
        hisin->Copy(*entryp[index]);

    }
}

void TCTWaveform::Draw(Int_t number,Option_t *option,Float_t low,Float_t high)
{
    //  Draw waveform
    //           Int_t number     ;  index of waveform
    //           Option_t *option ;  graphic option (see TH1F)
    // 	     Float_t low      ;  Time window
    // 	     Float_t high     ;
    if(number<0 || number>=Multiple) printf("No such measurment!!\n");
    else {
        if(low!=-1111. || high!=-1111.)
            ((TH1F *)histo->At(number))->GetXaxis()->SetRange(((TH1F *)histo->At(number))->GetXaxis()->FindBin(low),((TH1F *)histo->At(number))->GetXaxis()->FindBin(high));
        TH1 *his=((TH1F *)histo->At(number))->DrawCopy(option);

        his->SetLineColor(1);
        his->SetXTitle("t[ns]");
        his->SetYTitle("I [V/50#Omega]");
        his->SetLabelSize(0.045,"X");
        his->SetLabelSize(0.045,"Y");
        // Legend(number,number,number);
    }
}

void TCTWaveform::DrawMulti(Float_t low,Float_t high,Int_t Start,Int_t End,Int_t Step,Int_t model,Int_t indsub)
{
    //  Draws multiple waveform on the same plot with legend
    //           Float_t low      ; Time window
    //           Float_t high     ;
    // 	     Int_t Start      ; start index
    // 	     Int_t End        ; end index
    //	     Int_t Step       ; step (default=1)
    //         Int_t model      ; positive scale in volts=2, deconvoluted pulse=1
    //         Int_t indsub      ; substract wafeforw index
    // One of the most used functions
    Char_t v[5];
    Int_t i,maxi,mini,color,cii=0;
    Int_t colori[]={1,2,3,4,6,7,13,28,30,34,38,40,31,46,49,1,2,3,4,6,7,1,2,3,4,6,7,13,28,30,34,38,40,31,46,49,1,2,3,4,5,6,7,1,2,3,4,5,6,7,13,28,30,34,38,40,31,46,49,1,2,3,4,5,6,7};
    Float_t max=0,min=0,valmax,valmin,temp;

    TH1F *his,*hism,*subhis;
    if(Start==-1) Start=0;
    if(End==-1) End=Multiple-1;
    for(i=Start;i<=End;i+=Step) {
        valmax=((TH1F *)histo->At(i))->GetMaximum();
        valmin=((TH1F *)histo->At(i))->GetMinimum();
        if(valmax>max) {max=valmax; maxi=i;}
        if(valmin<min) {min=valmin; mini=i;}
    }
    if(max==0) maxi=End; //printf("maxi=%d",maxi);
    if(low!=-1111. || high!=-1111.)
        ((TH1F *)histo->At(maxi))->GetXaxis()->SetRange(((TH1F *)histo->At(maxi))->GetXaxis()->FindBin(low),((TH1F *)histo->At(maxi))->GetXaxis()->FindBin(high));

    hism=((TH1F *)histo->At(maxi));

    his=hism;


    his->SetXTitle("t[ns]");
    his->SetYTitle("I [V/50#Omega]");
    his->GetXaxis()->SetTitleOffset(1.1);
    his->GetYaxis()->SetTitleOffset(1.1);
    his->SetLabelSize(0.040,"X");
    his->SetLabelSize(0.040,"Y");
    int2ascii(v,Temperature.GetSum()/Temperature.GetSize());
    TString title="TCT Measurement @ T=";title=title+v; title=title+" C";
    his->SetTitle((const char *)title);
    if(TMath::Abs(max/min)>0.2) his->SetMinimum(min*1.1);
    if(DrawMode) his->DrawCopy(); else his->Draw();
    // for(int ss=2500;ss<3000;ss++) printf("%d %f\n",ss,his->GetBinContent(ss));
    Legend((TH1 *)his,Start,End,Step,model);

    if(indsub!=-1111)
    {
        subhis=((TH1F *)histo->At(indsub));
        subhis->Scale(-1.);
    }


    for(i=Start;i<=End;i+=Step)
    {
        //      color=i/7*40+i%7+1;
        hism=((TH1F *)histo->At(i));
        //     color=1;
        color=colori[cii];
        hism->SetLineColor((Color_t)color);
        his=hism;
        if(indsub!=-1111) {his->Add(subhis);}
        if(DrawMode) his->DrawCopy("SAME"); else his->Draw("SAME");
        cii++;
    }

}


Float_t TCTWaveform::Integral(Int_t i,Float_t mint, Float_t maxt)
{
    // Get waveform integral in time window
    // 	Int_t i      ; index
    //      Float_t mint ; Time window
    //      Float_t maxt ;
    // One of the most used functions
    Int_t mintime;
    Int_t maxtime;

    TH1F *his=((TH1F *)histo->At(i));
    if (maxt==-1111) maxtime=his->GetNbinsX()-1; else maxtime=his->GetXaxis()->FindBin(maxt);
    if (mint==-1111) mintime=1; else mintime=his->GetXaxis()->FindBin(mint);

    return(his->Integral(mintime,maxtime)*his->GetBinWidth(1));
}

Float_t TCTWaveform::Integral(Float_t y,Float_t mint, Float_t maxt)
{
    // Get waveform integral in time window (uses voltage insted of index)
    Int_t index;
    for(Int_t i=0;i<Multiple;i++) if(Voltages[i]==y) index=i;
    return(Integral(index,mint,maxt));
}



//___________________________________________________________________________________________________-
void TCTWaveform::GetIntegral(Float_t *x,Float_t *y,Float_t scale,Float_t mint, Float_t maxt)
{
    // Gets integrals in certain time window
    // 		Float_t  *x  ;  array of voltages
    // 		Float_t  *y  ;  array of integrals
    // 		Float_t  scale ;  scaling factor
    // 		Float_t mint;   Time window
    // 		Float_t maxt;

    Int_t i;
    Int_t mintime;
    Int_t maxtime;

    TH1F *his;
    for(i=0;i<Multiple;i++)
    {
        his=((TH1F *)histo->At(i));
        if (maxt==-1111) maxtime=his->GetNbinsX()-1; else maxtime=his->GetXaxis()->FindBin(maxt);
        if (mint==-1111) mintime=1; else mintime=his->GetXaxis()->FindBin(mint);
        //	printf("%d ,%d\n",mintime,maxtime);
        y[i]=his->Integral(mintime,maxtime)*his->GetBinWidth(1)/scale;
        if(Voltages[i]!=0 || Voltages[i+1]!=0) x[i]=Voltages[i]; else x[i]=(Float_t)i;
    }
}
//___________________________________________________________________________________________________-
void TCTWaveform::GetIntegral(Float_t *y,Float_t scale,Float_t mint, Float_t maxt)
{
    // Get Integral values only
    Float_t *x=new Float_t[Multiple];
    GetIntegral(x,y,scale,mint,maxt);
}
//___________________________________________________________________________________________________-




TGraph *TCTWaveform::CCE(Float_t mint, Float_t maxt,Int_t model,Int_t Show)
{
    // Get integral (charge) plot
    //                Float_t mint  ; Time Window
    //                Float_t maxt  ;
    //		    Int_t model   ; lin-lin=0  scale  sqrt-lin=1 abs(lin)-lin=2
    //		    Int_t Show    ; show graph
    Char_t v[6];
    Int_t i;
    Float_t *integral=new Float_t [Multiple];
    Float_t *index=new Float_t [Multiple];

    GetIntegral(index,integral,1,mint,maxt);

    if(model==1) for(i=0;i<Multiple;i++) index[i]=TMath::Sqrt(TMath::Abs(index[i]));
    if(model==2) for(i=0;i<Multiple;i++) index[i]=TMath::Abs(index[i]);

    TGraph *gr=new TGraph(Multiple,index,integral);
    gr->SetLineColor(2);
    gr->SetLineWidth(4);
    gr->SetMarkerColor(4);
    gr->SetMarkerStyle(21);

    //Get mean tempearature written in the graph title
    int2ascii(v,Temperature.GetSum()/Temperature.GetSize());
    TString title="Charge vs. Voltage @ T=";title=title+v; title=title+" C";
    gr->SetTitle((const char *)title);

    if(Show)
    {
        gr->Draw("APL");
        if(model!=1) gr->GetHistogram()->SetXTitle("U[V]"); else gr->GetHistogram()->SetXTitle("Sqrt U[ Sqrt V]");
        gr->GetHistogram()->SetYTitle("Charge[arb.]");
        gr->GetHistogram()->Draw();
        gr->Draw("APL");
    }

    delete integral;
    delete index;
    return(gr);
}

void TCTWaveform::Legend(TH1 *ch, Int_t start, Int_t end, Int_t Step,Int_t model)
{
    // Draw Legen (DrawMulti)
    //	       TH1F *ch    ; histogam
    //	       Int_t start ; start index
    //	       Int_t end   ; end index
    //	       Int_t model ; positive scale in volts=2
    Float_t minx,miny,maxy,maxx,x1,x2,y1,y2;
    Char_t title[30];

    Int_t color,cii=0;
    Int_t colori[]={1,2,3,4,6,7,13,28,30,34,38,40,31,46,49,1,2,3,4,6,7,1,2,3,4,6,7,13,28,30,34,38,40,31,46,49,1,2,3,4,5,6,7,1,2,3,4,5,6,7,13,28,30,34,38,40,31,46,49,1,2,3,4,5,6,7};
    TText *text;

    minx=ch->GetXaxis()->GetBinCenter(ch->GetXaxis()->GetFirst());
    maxx=ch->GetXaxis()->GetBinCenter(ch->GetXaxis()->GetLast());
    miny=ch->GetMinimum();
    maxy=ch->GetMaximum();

    x1=(maxx-minx)*0.6+minx;
    x2=(maxx-minx)*0.9+minx;
    y2=(maxy-miny)*0.35+miny;
    y1=(maxy-miny)*0.95+miny;

    //printf("coords: x1=%f y1=%f x2=%f y2=%f\n",x1,y1,x2,y2);

    if(pt!=NULL) delete pt;
    //pt=new TPaveText(x1,y1,x2,y2);
    for(Int_t i=start;i<=end;i+=Step)
    {
        if(model!=0)
            sprintf(title,"%s %d %s",prefix,(Int_t)TMath::Abs(GetVoltage(i)),suffix);
        else
            sprintf(title,"%s %d %s",prefix,(Int_t) GetVoltage(i),suffix);
        //   color=i/7*40+i%7+1;
        color=colori[cii];
        //       color=1;
        text=pt->AddText(title);
        text->SetTextColor(color);
        text->SetTextSize(0.05);
        cii++;
    }
    pt->Draw();
}

void TCTWaveform::Info()
{
    // shows information about measurements
    printf("Mesurement information:\n");
    printf("DATE: %d.%d.%d  TIME: %d:%d:%d\n",(Int_t)Date[0],(Int_t)Date[1],(Int_t)Date[2],(Int_t)Date[3],(Int_t)Date[4],(Int_t)Date[5]);
    printf("Frequency: %4.2f \n",Frequency);
    for(Int_t i=0;i<Multiple;i++) printf(" Temperature=%4.1f  Voltage=%4.1f  Current=%4.2e\n",Temperature[i],Voltages[i],Current[i]);
}

Float_t TCTWaveform::GetTime(Float_t scale)
{
    // Calculates the time of measurement
    //		Float_t scale ; scale (3600 = 1h)
    Float_t MonthSec[12]={0,31,59,90,120,151,181,212,243,273,303,333};
    Float_t Year=Date[2]*31536000./scale;
    if(((Int_t)Year)%4==0) for(Int_t i=2;i<12;i++) MonthSec[i]+=1;
    Float_t Month=MonthSec[(Int_t) (Date[1]-1)]*86400./scale;
    Float_t Day=Date[0]*86400./scale;;
    Float_t Hour=Date[3]*3600./scale;
    Float_t Min=Date[4]*60./scale;
    Float_t Sec=Date[5]/scale;
    //printf("%f %f %f %f %f %f\n",Year,Month,Day,Hour,Min,Sec);
    return(Year+Month+Day+Hour+Min+Sec);
}


void TCTWaveform::NormArray(Int_t num,Float_t *array)
{
    // array normalization
    Float_t max=0;
    Int_t maxi=0,i=0;
    for(i=0;i<num;i++) if(array[i]>max) {max=array[i]; maxi=i;}
    for(i=0;i<num;i++) array[i]/=max;
}

Int_t TCTWaveform::int2ascii(Char_t v[],Float_t vol,Int_t sign)
{
  // int to ascii with signs (same as in MeasureWF)
  Int_t k=0;
  if(sign) if(vol>0) v[k++]='+'; else {v[k++]='-'; vol=-vol;}
  if(!sign) vol=vol>0?vol:-vol;
  if((Int_t) vol/100!=0) v[k++]=(Char_t)(((Int_t) vol)/100)+48;
  v[k++]=(Char_t)(((Int_t) vol%100)/10)+48;
  v[k++]=(Char_t)((Int_t) vol%10)+48;
  v[k++]='\0';
  return k-1;
}

void TCTWaveform::SetHistoTime(int ind,Float_t *start)
{
    for(Int_t i=0;i<ind;i++)
        SetHistoTime(i,start[i]);

}

void TCTWaveform::SetHistoTime(int ind,Float_t start)
{
    TClonesArray &entry = *histo;
    TH1F *his=(TH1F *)entry[ind];
    Float_t minx=his->GetXaxis()->GetBinLowEdge(0)+start;
    Float_t maxx=his->GetXaxis()->GetBinUpEdge(his->GetNbinsX())+start;
    his->SetBins(his->GetNbinsX(),minx,maxx);
}
