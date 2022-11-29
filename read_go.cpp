#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <iterator>

using namespace std;

int main()
{

    string db,
        db_object_id,
        db_object_symbol,
        qualifier,
        go_id,
        db_reference,
        evidence_code,
        with_or_from,
        aspect,
        db_object_name,
        db_object_synonym,
        db_object_type,
        taxon,
        date,
        assigned_by,
        annotation_extension,
        gene_product_form_id;
    ifstream infile("/home/ad/gene_ontology/goa_uniprot_all.gaf");
    string line;

    ofstream outfile("test.tsv");

    while (getline(infile, line))
    {
        if (line[0] == '!')
            continue;
        istringstream iss(line);

        iss >> db >>
            db_object_id >>
            db_object_symbol >>
            qualifier >>
            go_id >>
            db_reference >>
            evidence_code >>
            with_or_from >>
            aspect >>
            db_object_name >>
            db_object_synonym >>
            db_object_type >>
            taxon >>
            date >>
            assigned_by >>
            annotation_extension >>
            gene_product_form_id;

        if ((db == "UniProtKB") && (qualifier == "enables") && (aspect == "F"))
            outfile << db_object_id << '\t' << go_id << endl;

    }
    outfile.close();
    return 0;
}
