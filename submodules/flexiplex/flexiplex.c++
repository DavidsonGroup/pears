// Copyright 2022 Nadia Davidson 
// This program is distributed under the MIT License.
// We also ask that you cite this software in publications
// where you made use of it for any part of the data analysis.

#include <iostream>
#include <fstream>
#include <istream>
#include <string>
#include <sstream>
#include <unordered_set>
#include <unordered_map>
#include <stdlib.h>
#include <unistd.h>
#include <vector>
#include <algorithm>
#include <thread>

#include "edlib.h"

#include <ctime>

using namespace std;

const static string VERSION="0.96.2";

// the help information which is printed when a user puts in the wrong
// combination of command line options.
void print_usage(){
  cerr << "usage: flexiplex [options] [reads_input]"  << endl;
  cerr << "  reads_input: a .fastq or .fasta file. Will read from stdin if empty." << endl;
  cerr << "  options: " << endl;
  cerr << "     -k known_list   Either 1) a text file of expected barcodes in the first column," << endl; 
  cerr << "                     one row per barcode, or 2) a comma separate string of barcodes. " << endl;
  cerr << "                     Without this option, flexiplex will search and report possible barcodes." << endl;
  cerr << "                     The generated list can be used for known_list in subsequent runs." << endl; 
  cerr << "     -i true/false   Replace read ID with barcodes+UMI, remove search strings" << endl;
  cerr << "                     including flanking sequenence and split read if multiple" << endl;
  cerr << "                     barcodes found (default: true)." << endl;
  cerr << "     -s true/false   Sort reads into separate files by barcode (default: false)" << endl;
  cerr << "     -l left     Left flank sequence to search for (default: CTACACGACGCTCTTCCGATCT)." << endl;
  cerr << "     -r right    Right flank sequence to search for (default: TTTTTTTTT)." << endl;
  cerr << "     -n prefix   Prefix for output filenames." << endl;
  cerr << "     -b N   Barcode length (default: 16)." << endl;
  cerr << "     -u N   UMI length (default: 12)." << endl;
  cerr << "     -e N   Maximum edit distance to barcode (default 2)." << endl;
  cerr << "     -f N   Maximum edit distance to primer+polyT (default 8)." << endl;
  cerr << "     -p N   Number of threads (default: 1)." << endl;
  cerr << "     -h     Print this usage information." << endl;
  cerr << endl;
}


// compliment nucleotides - used to reverse compliment string
char compliment(char& c){
  switch(c){
  case 'A' : return 'T';
  case 'T' : return 'A';
  case 'G' : return 'C';
  case 'C' : return 'G';
  default: return 'N';
  }
}

//Inplace reverse compliment
void reverse_compliment(string & seq){
   reverse(seq.begin(),seq.end());
   transform(seq.begin(),seq.end(),seq.begin(),compliment);
}

//Holds the search string patterns
struct SearchSeq {
  string primer;
  string polyA;
  string umi_seq;
  string temp_barcode;
} search_pattern ;

//Holds the found barcode and associated information 
struct Barcode {
  string barcode;
  string umi;
  int editd;
  int flank_editd;
  int flank_start;
  int flank_end;
} ;

struct SearchResult {
  string read_id;
  string qual_scores;
  string line;
  string rev_line;
  vector<Barcode> vec_bc_for;
  vector<Barcode> vec_bc_rev;
};

// Code for fast edit distance calculation for short sequences modified from
// https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance#C++
// s2 is always assumned to be the shorter string (barcode)
unsigned int edit_distance(const std::string& s1, const std::string& s2, unsigned int &end, int max_editd){
  
  std::size_t len1 = s1.size()+1, len2 = s2.size()+1;
  const char * s1_c = s1.c_str(); const char * s2_c = s2.c_str();

  vector< unsigned int> dist_holder(len1*len2);
  //initialise the edit distance matrix.
  //penalise for gaps at the start and end of the shorter sequence (j)
  //but not for shifting the start/end of the longer sequence (i,0)
  dist_holder[0]=0; //[0][0]
  for(unsigned int j = 1; j < len2; ++j) dist_holder[j] = j; //[0][j];
  for(unsigned int i = 1; i < len1; ++i) dist_holder[i*len2] = 0; //[i][0];
    
  int best=len2;
  end=len1-1;

  //loop over the distance matrix elements and calculate running distance
  for(unsigned int j = 1; j < len2; ++j){
    bool any_below_threshold=false; //flag used for early exit
    for(unsigned int i = 1; i < len1; ++i){
      int sub=(s1_c[i - 1] == s2_c[j - 1]) ? 0 : 1 ; //are the bases the same?
      if(sub==0) //if yes, no need to incremet distance
	dist_holder[i*len2+j]=dist_holder[(i-1)*len2+(j-1)];
      else //otherwise add insertion,deletion or substituation 
	dist_holder[i*len2+j] = std::min({ //j+i*len2  //set[i][j]
	    dist_holder[(i-1)*len2+j] + 1, //[i-1][j]
	    dist_holder[i*len2+(j-1)] + 1, //[i][j-1]
	    dist_holder[(i-1)*len2+(j-1)] + 1}); // ((s1_c[i - 1] == s2_c[j - 1]) ? 0 : 1) });
      if(dist_holder[i*len2+j] <= max_editd) any_below_threshold=true;
      if(j==(len2-1) && dist_holder[i*len2+j]<best){ //if this is the last row in j
	best=dist_holder[i*len2+j];                  //check if this is the best running score
	end=i;                                       //update the end position of alignment
      }
    }
    if(!any_below_threshold){ //early exit to save time.
      return(100);
    }
  }
  return best; //return edit distance
}


//Given a string 'seq' search for substring with primer and polyT sequence followed by
//a targeted search in the region for barcode
//Seaquence seearch is performed using edlib

Barcode get_barcode(string & seq,
		    unordered_set<string> *known_barcodes,
		    int flank_max_editd,
		    int barcode_max_editd){//,
  //		    SearchSeq & ss){
  
  const int OFFSET=5; //wiggle room in bases of the expected barcode start site to search.

  //initialise struct variables for return:
  Barcode barcode;
  barcode.editd=100; barcode.flank_editd=100;

  //initialise edlib configuration
  EdlibEqualityPair additionalEqualities[5] = {{'?','A'},{'?','C'},{'?','G'},{'?','T'},{'?','N'}};
  EdlibAlignConfig edlibConf = {flank_max_editd, EDLIB_MODE_HW, EDLIB_TASK_PATH, additionalEqualities, 5};

  //search for primer and ployT (barcode and umi as wildcards)
  string search_string=
    search_pattern.primer+
    search_pattern.temp_barcode+
    search_pattern.umi_seq+
    search_pattern.polyA;
  EdlibAlignResult result = edlibAlign(search_string.c_str(), search_string.length(), seq.c_str(), seq.length(), edlibConf);
  if(result.status != EDLIB_STATUS_OK | result.numLocations==0 ){
    edlibFreeAlignResult(result);
    return(barcode); // no match found - return
  } //fill in info about found primer and polyT location
  barcode.flank_editd=result.editDistance;
  barcode.flank_start=result.startLocations[0];
  barcode.flank_end=result.endLocations[0];
  edlibFreeAlignResult(result);
  
  //now we need to work out where the start of the barcode sequence is to extract it and match
  int primer_length_in_seq=search_pattern.primer.size(); //set the start of the barcode as primer length from start of primer
  if(search_pattern.primer!=""){ //search for the primer to refine (deal with indel errors)
    string search_region=seq.substr(barcode.flank_start,search_pattern.primer.length()+OFFSET);
    EdlibAlignResult result_targeted = edlibAlign(
						  search_pattern.primer.c_str(),
						  search_pattern.primer.length(),
						  search_region.c_str(),
						  search_region.length(), edlibConf);
    if(result_targeted.status == EDLIB_STATUS_OK && result_targeted.numLocations > 0 ){
      primer_length_in_seq=result_targeted.endLocations[0]+1; //new barcode starting point
    }
    edlibFreeAlignResult(result_targeted);
  }
  
  //if not checking against known list of barcodes, return sequence after the primer
  //also check for a perfect match straight up as this will save computer later.
  int BC_start=barcode.flank_start+primer_length_in_seq;
  string exact_bc=seq.substr(BC_start,search_pattern.temp_barcode.length());
  if(known_barcodes->size()==0 || (known_barcodes->find(exact_bc) != known_barcodes->end())){ 
    barcode.barcode=exact_bc;
    barcode.editd=0;
    barcode.umi=seq.substr(BC_start+search_pattern.temp_barcode.length(),search_pattern.umi_seq.length());//resul
    return(barcode);
  }
  
  // otherwise widen our search space and the look for matches with errors
  string barcode_seq=seq.substr(BC_start-OFFSET,search_pattern.temp_barcode.length()+2*OFFSET);
  
  //iterate over all the known barcodes, checking each sequentially
  unordered_set<string>::iterator known_barcodes_itr=known_barcodes->begin();
  unsigned int editDistance; unsigned int endDistance;
  for(; known_barcodes_itr!=known_barcodes->end(); known_barcodes_itr++){
    search_string=(*known_barcodes_itr); //known barcode to check again
    editDistance = edit_distance(barcode_seq,search_string,endDistance,barcode_max_editd);
    if(editDistance < barcode.editd && editDistance <= barcode_max_editd){ //if best so far, update
      barcode.editd=editDistance; 
      barcode.barcode=*known_barcodes_itr;
      barcode.umi=seq.substr(BC_start-OFFSET+endDistance,search_pattern.umi_seq.length());//assumes no error in UMI seq.
      if(editDistance==0){ //if perfect match is found we're done.
	return(barcode);
      }
    }
  }
  return(barcode); //return the best matched barcode and associated information
}

//search a read for one or more barcodes (parent function that calls get_barcode)
vector<Barcode> big_barcode_search(string & sequence, unordered_set<string> & known_barcodes,
				   int max_flank_editd, int max_editd){ //, SearchSeq & ss){
  vector<Barcode> return_vec; //vector of all the barcodes found

  //search for barcode
  Barcode result=get_barcode(sequence,&known_barcodes,max_flank_editd,max_editd); //,ss);
  if(result.editd<=max_editd) //add to return vector if edit distance small enough
    return_vec.push_back(result);
  
  //if a result was found, mask out the flanking sequence and search again in case there are more.
  if(return_vec.size()>0){
    string masked_sequence = sequence;
    for(int i=0; i<return_vec.size(); i++){
      int flank_length=return_vec.at(i).flank_end-return_vec.at(i).flank_start;
      masked_sequence.replace(return_vec.at(i).flank_start,flank_length,string(flank_length,'X'));
    } //recursively call this function until no more barcodes are found
    vector<Barcode> masked_res;
    masked_res=big_barcode_search(masked_sequence,known_barcodes,max_flank_editd,max_editd); //,ss);
    return_vec.insert(return_vec.end(),masked_res.begin(),masked_res.end()); //add to result
  }
    
  return(return_vec);
}

// utility function to check true/false input options
bool get_bool_opt_arg(string value){
  transform(value.begin(), value.end(), value.begin(), ::tolower);
  if( value.compare("true")==0 | value.compare("t")==0 | value.compare("1")==0){
    return true;
  } else if (value.compare("false")!=0 | value.compare("f")!=0 | value.compare("0")!=0){
    return false;
  } else {
    cerr << "Unknown argument to boolean option" << endl;
    print_usage();
    exit(1);
  } 
}

// print information about barcodes
void print_stats(string read_id, vector<Barcode> & vec_bc, ostream & out_stream){
  for(int b=0; b<vec_bc.size() ; b++){
    out_stream << read_id << "\t"
	       << vec_bc.at(b).barcode << "\t"
	       << vec_bc.at(b).flank_editd << "\t"
	       << vec_bc.at(b).editd << "\t"
	       << vec_bc.at(b).umi << "\t"
	       << endl;
  }
}

void print_line(string id, string read, string quals, ostream & out_stream){

  //flag for read format
  bool is_fastq=!(quals==""); //no quality scores passed = fasta

  //output to the new read file
    if(is_fastq)
      out_stream << "@" << id << endl;
    else
      out_stream << ">" << id << endl;
    out_stream << read << endl;
    if(is_fastq){
      out_stream << "+" << id << endl;
      out_stream << quals << endl;
    }
}

//print fastq or fasta lines..
void print_read(string read_id,string read, string qual,
		vector<Barcode> & vec_bc, string prefix,
		bool split, unordered_set<string> & found_barcodes,
		bool trim_barcodes){
  //loop over the barcodes found... usually will just be one
  for(int b=0; b<vec_bc.size() ; b++){
    
    //format the new read id. Using FLAMES format.
    stringstream ss;
    ss << (b+1) << "of" << vec_bc.size() ;
    string barcode=vec_bc.at(b).barcode;
    string new_read_id=barcode+"_"+vec_bc.at(b).umi+"#"+read_id+ss.str();
    
    // work out the start and end base in case multiple barcodes
    int read_start=vec_bc.at(b).flank_end;
    int read_length=read.length()-read_start;
    for(int f=0; f<vec_bc.size(); f++){
      int temp_read_length=vec_bc.at(f).flank_start-read_start;
      if(temp_read_length>0 && temp_read_length<read_length)
	read_length=temp_read_length;
    }
    string qual_new=""; //don't trim the quality scores if it's a fasta file
    if(qual!="") qual_new=qual.substr(read_start,read_length);
    string read_new=read.substr(read_start,read_length);

    if(b==0 && !trim_barcodes){ //override if read shouldn't be cut
      new_read_id=read_id;
      read_new=read;
      qual_new=qual;
      b=vec_bc.size(); //force loop to exit after this iteration
    }
    
    if(split){ //to a file if spliting by barcode
      string outname=prefix+"_"+barcode+".";
      if(qual=="") outname+="fasta"; else outname+="fastq";
      ofstream outstream;
      if(found_barcodes.insert(barcode).second)
	outstream.open(outname); //override file if this is the first read for the barcode
      else
	outstream.open(outname,ofstream::app); //remove file if this is the first read for the barcode
      print_line(new_read_id,read_new,qual_new,outstream);
      outstream.close();
    } else {
      print_line(new_read_id,read_new,qual_new,std::cout);
    }
  }
}

// separated out from main so that this can be run with threads
void search_read(vector<SearchResult> & reads, unordered_set<string> & known_barcodes,
			 int flank_edit_distance, int edit_distance){
  
  for(int r=0; r<reads.size(); r++){
    
    //forward search
    reads[r].vec_bc_for=big_barcode_search(reads[r].line,
					   known_barcodes,
					   flank_edit_distance,
					   edit_distance);
    reads[r].rev_line=reads[r].line;
    reverse_compliment(reads[r].rev_line);
    //Check the reverse compliment of the read
    reads[r].vec_bc_rev=big_barcode_search(reads[r].rev_line,
					   known_barcodes,
					   flank_edit_distance,
					   edit_distance);
  }
}

// MAIN is here!!
int main(int argc, char **argv){
  std::ios_base::sync_with_stdio(false);

  cerr << "FLEXIPLEX " << VERSION << endl;

  //Variables to store user options
  //Set these to their defaults
  int expected_cells=0; //(d)
  int edit_distance=2; //(e)
  int flank_edit_distance=8; //(f)

  //set the output filenames
  string out_stat_filename="reads_barcodes.txt";
  string out_bc_filename="barcodes_counts.txt";
  string out_filename_prefix="flexiplex"; //(n)

  bool split_file_by_barcode=false; //(s)
  bool remove_barcodes=true; //(r)
  
  search_pattern.primer = "CTACACGACGCTCTTCCGATCT"; //(p)
  search_pattern.polyA = string(9,'T'); //(T)
  search_pattern.umi_seq = string(12,'?'); //(length u)
  search_pattern.temp_barcode = string(16,'?'); //(length b)
  
  //Set of known barcodes 
  unordered_set<string> known_barcodes;
  unordered_set<string> found_barcodes;

  // threads
  int n_threads=1;
  
  /*** Pass command line option *****/
  int c;
  int params=1;
  ifstream file;
  string line;

  while((c =  getopt(argc, argv, "k:i:l:r:b:u:e:f:n:s:h:p:")) != EOF){
    switch(c){
    case 'k': { //k=list of known barcodes
      string file_name(optarg);
      string bc;
      /**** READ BARCODE LIST FROM FILE ******/
      file.open(file_name);
      cerr << "Setting known barcodes from "<< file_name << endl;
      if(!(file.good())){ //if the string given isn't a file
	stringstream bc_list(file_name); string s;
	while (getline(bc_list, bc, ',')) //tokenize
	  known_barcodes.insert(bc);
      } else {
	// otherwise get the barcodes from the file..
	while ( getline (file,line) ){
	  istringstream line_stream(line);
	  line_stream >> bc;
	  known_barcodes.insert(bc); 
	}
	file.close();
      }
      cerr << "Number of known barcodes: " << known_barcodes.size() << endl;
      if(known_barcodes.size()==0){
	print_usage();
	exit(1); //case barcode file is empty
      }
      //set barcode length automatically from known barcodes..
      int bl=(known_barcodes.begin())->length();
      search_pattern.temp_barcode=string(bl,'?');
      cerr << "Setting barcode length automatically to " << bl << endl;
      params+=2;
      break;     
    }
    case 'i':{
      remove_barcodes=get_bool_opt_arg(optarg);
      cerr << "Setting read IDs to be replaced: "<< remove_barcodes << endl;
      params+=2;
      break;
    }
    case 'e':{
      edit_distance=atoi(optarg);
      cerr << "Setting max barcode edit distance to "<< edit_distance << endl;
      params+=2;
      break;
    }
    case 'f':{
      flank_edit_distance=atoi(optarg);
      cerr << "Setting max flanking sequence edit distance to "<< flank_edit_distance << endl;
      params+=2;
      break;
    }
    case 'l':{
      search_pattern.primer=optarg;
      cerr << "Setting primer to search for: " << search_pattern.primer << endl;
      params+=2;
      break;
    }
    case 'r':{
      search_pattern.polyA=optarg;
      cerr << "Setting polyT to search for: " << search_pattern.polyA << endl;
      params+=2;
      break;
    }
    case 'u':{
      int ul=atoi(optarg);
      search_pattern.umi_seq=string(ul,'?');
      cerr << "Setting UMI length to " << ul << endl;
      params+=2;
      break;
    }
    case 'b':{
      int bl=atoi(optarg);
      search_pattern.temp_barcode=string(bl,'?');
      cerr << "Setting barcode length to " << bl << endl;
      params+=2;
      break;
    }
    case 'h':{
      print_usage();
      exit(1);
    }
    case 'n':{
      out_filename_prefix=optarg;
      cerr << "Setting output filename prefix to: " << out_filename_prefix << endl;
      params+=2;
      break;
    }
    case 's':{
      split_file_by_barcode=get_bool_opt_arg(optarg);
      cerr << "Split read output into separate files by barcode: " << split_file_by_barcode << endl;
      int max_split_bc=50;
      if(known_barcodes.size()>max_split_bc){
	cerr << "Too many barcodes to split into separate files: "<< known_barcodes.size()
	     << "> "<< max_split_bc<< endl;
	split_file_by_barcode=false;
      }
      params+=2;
      break;
    }
    case 'p':{
      n_threads=atoi(optarg);
      cerr << "Setting number of threads to "<< n_threads << endl;
      params+=2;
      break;
    }
    case '?': //other unknown options
      cerr << "Unknown option.. stopping" << endl;
      print_usage();
      exit(1);
    }
  }
  
  cerr << "For usage information type: flexiplex -h" << endl;
  
  istream * in;
  FILE * ifile;
    
  //check that a read file is given
  if(params>=argc){
    cerr << "No filename given... getting reads from stdin..." << endl;
    in=&std::cin;
  } else {
    // check that the reads fileis okay
    string reads_file=argv[params];
    file.open(reads_file);
    if(!(file.good())){
      cerr << "Unable to open file " << reads_file << endl;
      print_usage();
      exit(1);
    }
    in=&file;
  }
  
  /********* FIND BARCODE IN READS ********/
  string sequence;
  int bc_count=0;
  int r_count=0;
  int multi_bc_count=0;
  
  ofstream out_stat_file;
  out_stat_filename=out_filename_prefix+"_"+out_stat_filename;
  out_bc_filename=out_filename_prefix+"_"+out_bc_filename;
  params+=2;

  if(known_barcodes.size()>0){
    out_stat_file.open(out_stat_filename);
    out_stat_file << "Read\tCellBarcode\tFlankEditDist\tBarcodeEditDist\tUMI"<<endl;
  }
  cerr << "Searching for barcodes..." << endl;
  bool is_fastq=true;
  unordered_map< string, int > barcode_counts; 
  string read_id_line;
  if(getline (*in,read_id_line)){ //check the first line for file type
    if(read_id_line[0]=='>'){ is_fastq=false;
    } else if (read_id_line[0] == '@'){ //fasta
    } else {
      cerr << "Unknown read format... exiting" << endl; exit(1);
    }
  }
  
  while ( getline (*in,line) ){
    const int buffer_size = 2000; //number of reads to pass to one thread.
    vector<vector<SearchResult>> sr_v(n_threads);
    for(int i=0; i<n_threads; i++)
      sr_v[i]=vector<SearchResult>(buffer_size); 
    vector<thread> threads(n_threads);
    for(int t=0; t < n_threads; t++){ //get n_threads*buffer number or reads..
      for(int b=0 ; b < buffer_size ; b++){ 
	SearchResult & sr = sr_v[t][b];
	sr.line=line;
	string read_id;
	//sr.read_id= read_id_line.substr(1,read_id_line.find_first_not_of(" \t")-1);      
	istringstream line_stream(read_id_line);
	line_stream >> sr.read_id;
	sr.read_id.erase(0,1);

	    
	//      string qual_scores="";
	if(!is_fastq){ //fastq (account for multi-lines per read)
	  string buffer_string; 
	  while(getline(*in,buffer_string) && buffer_string[0]!='>')
	    sr.line+=buffer_string;
	  read_id_line=buffer_string;
	} else { //fastq (get quality scores)
	  for(int s=0; s<2; s++) getline(*in,sr.qual_scores);
	  getline(*in,read_id_line);
	}
      
	r_count++; //progress counter
	if(r_count % 100000 == 0)
	  cerr << r_count/((double) 1000000 ) << " million reads processed.." << endl;

	//sr.read_id=read_id;
	//sr.line=line;
	//sr.qual_scores=qual_scores;
	
	//      sr_v[t][b] / buffer_size].push_back(sr);
      
	//this is quite ugly, must be a better way to do this..
	if(b==buffer_size-1 && t==n_threads-1){
	  break; //if it's the last in the chunk don't getline as this happens in the while statement
	} else if( !getline(*in,line)){ //case we are at the end of the reads.
	/**	std::time_t t_ = std::time(0);   // get time now
	std::tm* now = std::localtime(&t_);
	cout << "Read="<<t<< endl;
	cout << "2. Added read for thread" << (t / buffer) << " at " << asctime(now) << endl; **/
	  sr_v[t].resize(b+1);
	  threads[t]=std::thread(search_read,ref(sr_v[t]),ref(known_barcodes),flank_edit_distance,edit_distance);
	  for(int t2=t+1; t2 < n_threads ; t2++) sr_v[t2].resize(0);
	  goto print_result; //advance the line
	}
      }
      // send reads to the thread
      threads[t]=std::thread(search_read,ref(sr_v[t]),ref(known_barcodes),flank_edit_distance,edit_distance);
    }
    /**    t_ = std::time(0);   // get time now
    now = std::localtime(&t_);
    cout << "START-3" << asctime(now) << endl; 
    **/
  print_result:
    
    for(int t=0; t < sr_v.size(); t++){ //loop over the threads and print out ther results
      if(sr_v[t].size()>0) threads[t].join(); // wait for the threads to finish before printing
      
      for(int r=0; r< sr_v[t].size(); r++){ // loop over the reads
	
	for(int b=0; b<sr_v[t][r].vec_bc_for.size(); b++)
	  barcode_counts[sr_v[t][r].vec_bc_for.at(b).barcode]++;
	for(int b=0; b<sr_v[t][r].vec_bc_rev.size(); b++)
	  barcode_counts[sr_v[t][r].vec_bc_rev.at(b).barcode]++;
	
	if((sr_v[t][r].vec_bc_for.size()+sr_v[t][r].vec_bc_rev.size())>0)
	  bc_count++;
	if((sr_v[t][r].vec_bc_for.size()+sr_v[t][r].vec_bc_rev.size())>1 ){
	  multi_bc_count++;
	}
	
	if(known_barcodes.size()!=0){ // if we are just looking for all possible barcodes don't output reads etc.
	  
	  print_stats(sr_v[t][r].read_id, sr_v[t][r].vec_bc_for, out_stat_file);
	  print_stats(sr_v[t][r].read_id, sr_v[t][r].vec_bc_rev, out_stat_file);
	  
	  print_read(sr_v[t][r].read_id+"_+",sr_v[t][r].line,sr_v[t][r].qual_scores,sr_v[t][r].vec_bc_for,
		     out_filename_prefix,split_file_by_barcode,found_barcodes,remove_barcodes);
	  reverse(sr_v[t][r].qual_scores.begin(),sr_v[t][r].qual_scores.end());
	  if(remove_barcodes || sr_v[t][r].vec_bc_for.size()==0) //case we just want to print read once if multiple bc found.
	    print_read(sr_v[t][r].read_id+"_-",sr_v[t][r].rev_line,sr_v[t][r].qual_scores,sr_v[t][r].vec_bc_rev,
		       out_filename_prefix,split_file_by_barcode,found_barcodes,remove_barcodes);
	}
      }
    }
  }
  file.close();
  
  cerr << "Number of reads processed: " << r_count << endl;
  cerr << "Number of reads where a barcode was found: " << bc_count << endl;
  cerr << "Number of reads where more than one barcode was found: " << multi_bc_count << endl;
  cerr << "All done!" << endl;

  if(known_barcodes.size()>0){
    out_stat_file.close();
    return(0);
  }
  
  if(barcode_counts.size()==0)
    return(0);
  
  typedef std::pair<std::string, int> pair;
  vector<pair> bc_vec;
  copy(barcode_counts.begin(),barcode_counts.end(), back_inserter<vector<pair>>(bc_vec));
  sort(bc_vec.begin(), bc_vec.end(),[](const pair &l, const pair &r){
    if (l.second != r.second)
      return l.second > r.second;
    return l.first < r.first;
  });
  vector<int> hist(bc_vec[0].second);
  ofstream out_bc_file;
  out_bc_file.open(out_bc_filename);
  for (auto const &bc_pair: bc_vec){
    out_bc_file << bc_pair.first << "\t" << bc_pair.second << endl;
    hist[bc_pair.second-1]++;
  }
  out_bc_file.close();

  cout << "Reads\tBarcodes" << endl;
  for(int i=hist.size()-1; i>=0; i--)
    cout << i+1 << "\t" << hist[i] << endl;
    
}
