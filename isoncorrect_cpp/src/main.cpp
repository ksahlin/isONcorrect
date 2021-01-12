#include <iostream>
#include "argparse.hpp"
#include "edlib/edlib.h"
#include "spoa/spoa.hpp"


void run_edlib(const std::string x, const std::string y) {
  const char * ptr_x = x.c_str();
  const char * ptr_y = y.c_str();
  EdlibAlignResult result = edlibAlign(ptr_x, x.length(), ptr_y, y.length(), edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
  if (result.status == EDLIB_STATUS_OK) {
      printf("edit_distance('hello', 'world!') = %d\n", result.editDistance);
      printf("edit_distance('hello', 'world!') = %s\n", result.alignment);
      char* cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
      printf("%s\n", cigar);
      free(cigar);
  }
  edlibFreeAlignResult(result);
  std::cout << "edlib just got executed!\n";
}

void run_spoa() {
  std::vector<std::string> sequences = {
      "CATAAAAGAACGTAGGTCGCCCGTCCGTAACCTGTCGGATCACCGGAAAGGACCCGTAAAGTGATAATGAT",
      "ATAAAGGCAGTCGCTCTGTAAGCTGTCGATTCACCGGAAAGATGGCGTTACCACGTAAAGTGATAATGATTAT",
      "ATCAAAGAACGTGTAGCCTGTCCGTAATCTAGCGCATTTCACACGAGACCCGCGTAATGGG",
      "CGTAAATAGGTAATGATTATCATTACATATCACAACTAGGGCCGTATTAATCATGATATCATCA",
      "GTCGCTAGAGGCATCGTGAGTCGCTTCCGTACCGCAAGGATGACGAGTCACTTAAAGTGATAAT",
      "CCGTAACCTTCATCGGATCACCGGAAAGGACCCGTAAATAGACCTGATTATCATCTACAT"
  };

  auto alignment_engine = spoa::AlignmentEngine::Create(
      spoa::AlignmentType::kNW, 3, -5, -3);  // linear gaps

  spoa::Graph graph{};

  for (const auto& it : sequences) {
    auto alignment = alignment_engine->Align(it, graph);
    graph.AddAlignment(alignment, it);
  }

  auto consensus = graph.GenerateConsensus();

  std::cerr << ">Consensus LN:i:" << consensus.size() << std::endl
            << consensus << std::endl;

  auto msa = graph.GenerateMultipleSequenceAlignment();

  for (const auto& it : msa) {
    std::cerr << it << std::endl;
  }
  std::cout << "spoa just got executed!\n";
}

int main(int argc, char *argv[]) {
  argparse::ArgumentParser program("isONcorrect", "0.0.6");

  // Optional arguments (for backwards compatability with the python version)
  program.add_argument("--fastq")
    .help("Kmer size")
    .action([](const std::string& value) { return std::stoi(value); });

  // program.add_argument("k")
  //   .help("Kmer size")
  //   .action([](const std::string& value) { return std::stoi(value); });

  program.add_argument("-o", "--output")
    .default_value(std::string("."))
    .required()
    .help("specify the output folder destination.");

  program.add_argument("--k")
    .help("Kmer size")
    .default_value(int(9))
    .action([](const std::string& value) { return std::stoi(value); });

  program.add_argument("--w")
    .help("Window size")
    .default_value(int(10))
    .action([](const std::string& value) { return std::stoi(value); });;

  program.add_argument("--t")
    .help("Number of cores allocated for correction")
    .default_value(int(8))
    .action([](const std::string& value) { return std::stoi(value); });;
  

  // flags
  program.add_argument("--exact")
    .help("Get exact solution for WIS for evary read \
    (recalculating weights for each read (much slower but slightly more accuracy,\
    not to be used for clusters with over ~500 reads)")
    .default_value(false)
    .implicit_value(true);

  try {
    program.parse_args(argc, argv);
  }
  catch (const std::runtime_error& err) {
    std::cout << err.what() << std::endl;
    std::cout << program;
    exit(0);
  }
  
  auto k = program.get<int>("--k");
  std::cout << (k * k) << std::endl;
  auto w = program.get<int>("--w");
  std::cout << "w is set to " << w << std::endl;
  std::string x = "hello";
  std::string y = "world!";
  run_edlib(x, y);

  run_spoa();

  return 0;
}

    // parser.add_argument('--version', action='version', version='%(prog)s 0.0.6')

    // parser.add_argument('--fastq', type=str,  default=False, help='Path to input fastq file with reads')
    // # parser.add_argument('--t', dest="nr_cores", type=int, default=8, help='Number of cores allocated for clustering')

    // parser.add_argument('--k', type=int, default=9, help='Kmer size')
    // parser.add_argument('--w', type=int, default=10, help='Window size')
    // parser.add_argument('--xmin', type=int, default=14, help='Upper interval length')
    // parser.add_argument('--xmax', type=int, default=80, help='Lower interval length')
    // parser.add_argument('--T', type=float, default=0.1, help='Minimum fraction keeping substitution')
    // # parser.add_argument('--C', type=float, default=0.05, help='Minimum fraction of keeping alternative refernece contexts')
    // parser.add_argument('--exact', action="store_true", help='Get exact solution for WIS for evary read (recalculating weights for each read (much slower but slightly more accuracy,\
    //                                                              not to be used for clusters with over ~500 reads)')
    // parser.add_argument('--disable_numpy', action="store_true", help='Do not require numpy to be installed, but this version is about 1.5x slower than with numpy.')

    // parser.add_argument('--max_seqs_to_spoa', type=int, default=200,  help='Maximum number of seqs to spoa')
    // parser.add_argument('--max_seqs', type=int, default=1000,  help='Maximum number of seqs to correct at a time (in case of large clusters).')
    // parser.add_argument('--use_racon', action="store_true", help='Use racon to polish consensus after spoa (more time consuming but higher accuracy).')

    // parser.add_argument('--exact_instance_limit', type=int, default=0,  help='Activates slower exact mode for instance smaller than this limit')
    // # parser.add_argument('--w_equal_k_limit', type=int, default=0,  help='Sets w=k which is slower and more memory consuming but more accurate and useful for smalled clusters.')
    // parser.add_argument('--set_w_dynamically', action="store_true", help='Set w = k + max(2*k, floor(cluster_size/1000)).')
    // parser.add_argument('--verbose', action="store_true", help='Print various developer stats.')

    // parser.add_argument('--compression', action="store_true", help='Use homopolymenr compressed reads. (Deprecated, because we will have fewer \
    //                                                                     minmimizer combinations to span regions in homopolymenr dense regions. Solution \
    //                                                                     could be to adjust upper interval legnth dynamically to guarantee a certain number of spanning intervals.')
    // parser.add_argument('--outfolder', type=str,  default=None, help='A fasta file with transcripts that are shared between samples and have perfect illumina support.')
    // # parser.add_argument('--pickled_subreads', type=str, help='Path to an already parsed subreads file in pickle format')
    // # parser.set_defaults(which='main')
    // args = parser.parse_args()


    // if args.xmin < 2*args.k:
    //     args.xmin = 2*args.k
    //     eprint("xmin set to {0}".format(args.xmin))

    // if len(sys.argv)==1:
    //     parser.print_help()
    //     sys.exit()
    // if not args.fastq and not args.flnc and not  args.ccs:
    //     parser.print_help()
    //     sys.exit()




    // if args.outfolder and not os.path.exists(args.outfolder):
    //     os.makedirs(args.outfolder)


    // # edlib_module = 'edlib'
    // # parasail_module = 'parasail'
    // # if edlib_module not in sys.modules:
    // #     print('You have not imported the {0} module. Only performing clustering with mapping, i.e., no alignment.'.format(edlib_module))
    // # if parasail_module not in sys.modules:
    // #     eprint('You have not imported the {0} module. Only performing clustering with mapping, i.e., no alignment!'.format(parasail_module))
    // #     sys.exit(1)
    // if 100 < args.w or args.w < args.k:
    //     eprint('Please specify a window of size larger or equal to k, and smaller than 100.')
    //     sys.exit(1)
