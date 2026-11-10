// Stream SAM from samtools, audit every record, and feed two samtools checksums.
// No BAM parsing, record editing, or custom checksum algorithm is implemented.
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <stdexcept>
#include <unistd.h>
using namespace std;
vector<string> fields(const string& s) {
    vector<string> out; size_t a=0,b;
    while ((b=s.find('\t',a))!=string::npos) {out.push_back(s.substr(a,b-a));a=b+1;}
    out.push_back(s.substr(a));return out;
}
bool cigar(const string& s) {
    if(s.empty()||s=="*")return false;
    long n=0; bool digits=false;
    for(char c:s) {if(c>='0'&&c<='9'){n=n*10+c-'0';digits=true;}
        else {if(!digits||n<=0||string("MIDNSHP=X").find(c)==string::npos)return false;n=0;digits=false;}}
    return !digits;
}
int main(int argc,char** argv) {
  try {
    if(argc!=5)throw runtime_error("audit <sample> <json> <mapped-fd> <coordinate:0|1>");
    string sample=argv[1],line; FILE* mapped=fdopen(stoi(argv[3]),"w");
    if(!mapped)throw runtime_error("mapped checksum pipe unavailable");
    setvbuf(stdout,nullptr,_IOFBF,1048576);setvbuf(mapped,nullptr,_IOFBF,1048576);
    map<string,string> groups;map<string,int> tids;map<string,long long> n;
    for(string k:{"records","primary","secondary","supplementary","mc_required","missing_mc","invalid_mc",
                  "missing_or_invalid_rg","hard_clipped_supplementary","missing_primary_sequence",
                  "missing_primary_quality","coordinate_errors"})n[k]=0;
    long long lastpos=-1;int lasttid=-1;
    while(getline(cin,line)) {
      auto f=fields(line);
      if(line.rfind("@SQ\t",0)==0){for(auto v:f)if(v.rfind("SN:",0)==0)tids[v.substr(3)]=tids.size();}
      if(line.rfind("@RG\t",0)==0){string id,sm;for(auto v:f){if(v.rfind("ID:",0)==0)id=v.substr(3);if(v.rfind("SM:",0)==0)sm=v.substr(3);}groups[id]=sm;}
      bool header=!line.empty()&&line[0]=='@';
      if(!header){
        if(f.size()<11)throw runtime_error("Malformed SAM record");
        int flag=stoi(f[1]);n["records"]++;
        bool secondary=flag&256,supp=flag&2048;
        if(secondary)n["secondary"]++;else if(supp)n["supplementary"]++;else n["primary"]++;
        string rg,mc;
        for(size_t i=11;i<f.size();i++){if(f[i].rfind("RG:Z:",0)==0)rg=f[i].substr(5);if(f[i].rfind("MC:Z:",0)==0)mc=f[i].substr(5);}
        if(!groups.count(rg)||groups[rg]!=sample)n["missing_or_invalid_rg"]++;
        if(!secondary&&(flag&1)&&!(flag&4)&&!(flag&8)){
          n["mc_required"]++;if(mc.empty())n["missing_mc"]++;else if(!cigar(mc))n["invalid_mc"]++;
        }
        if(supp&&f[5].find('H')!=string::npos)n["hard_clipped_supplementary"]++;
        if(!secondary&&!supp){if(f[9]=="*")n["missing_primary_sequence"]++;if(!(flag&4)&&f[10]=="*")n["missing_primary_quality"]++;}
        int tid=f[2]=="*"?tids.size():tids.count(f[2])?tids[f[2]]:-1;
        long long pos=stoll(f[3]);
        if(tid<0||tid<lasttid||(tid==lasttid&&pos<lastpos))n["coordinate_errors"]++;
        lasttid=tid;lastpos=pos;
        if(!(flag&4)){if(fwrite(line.data(),1,line.size(),mapped)!=line.size()||fputc('\n',mapped)==EOF)throw runtime_error("mapped pipe failed");}
      }else{if(fwrite(line.data(),1,line.size(),mapped)!=line.size()||fputc('\n',mapped)==EOF)throw runtime_error("mapped header pipe failed");}
      if(fwrite(line.data(),1,line.size(),stdout)!=line.size()||fputc('\n',stdout)==EOF)throw runtime_error("checksum pipe failed");
    }
    if(cin.bad()||fclose(mapped)||fflush(stdout))throw runtime_error("stream read/write error");
    ofstream out(argv[2]);out<<"{\n";bool first=true;for(auto p:n){if(!first)out<<",\n";first=false;out<<"  \""<<p.first<<"\": "<<p.second;}out<<"\n}\n";
    if(!out)throw runtime_error("audit JSON write failed");
    return 0;
  }catch(exception& e){cerr<<e.what()<<endl;return 1;}
}
