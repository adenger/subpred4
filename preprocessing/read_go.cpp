#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <iterator>

using namespace std;

// Requires ~300GB of RAM
int main()
{
    ifstream infile("/home/ad/gene_ontology/goa_uniprot_all.gaf");
    string line, element;

    ofstream outfile("goa_sp_iea_mf_enables2.tsv");

    cout << "reading file..." << endl;
    vector<string> *lines = new vector<string>();
    while (getline(infile, line))
    {
        lines->push_back(line);
    }
    cout << "writing output..." << endl;
    for (string line : *lines)
    {
        if (line[0] == '!')
        {
            continue;
        }
        vector<string> elements;
        istringstream line_stream(line);
        while (getline(line_stream, element, '\t'))
        {
            elements.push_back(element);
        }
       
        if (elements[3] == "enables" && elements[8] == "F" && elements[0] == "UniProtKB")
        {
            outfile << elements[1] << '\t' << elements[4] << '\t' << elements[6] << endl;
        }
    }
    delete lines;
    return 0;
}
