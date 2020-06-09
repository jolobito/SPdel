import { csvParseRows } from 'd3-dsv';

function DataParser() {
    this.data = {
        individuals: [],
        distances: []
    };

}

DataParser.prototype = {
    parseName: function(name){
        const individual = {};
        const split = name.split('_');
        individual.id = name;
        if(split[0] === 'MOTU'){
            individual.groupping = split[0];
            individual.group =  split[0] + " " + split[1];
            individual.genre = split[2];
            individual.species = split[3];
            individual.code = split[4];
        }else{
            individual.groupping = split[0];
            individual.group =  split[0] + " " + split[1];
            individual.genre = split[0];
            individual.species = split[1];
            individual.code = split[3];
        }
        return individual;
    },
    parseHeader: function(row, data){
        for (let j = 1; j < row.length; ++j) {
            this.data.individuals.push(this.parseName(row[j]));
        }
        return data;
    },
    parseDataRow : function(row){
        const id_i = row[0];
        for (let j = 1; j < row.length; ++j) {
            const id_j = this.data.individuals[j-1].id;
            const d_ij = {
                i: id_i,
                j: id_j,
                distance: row[j]
            };
            this.data.distances.push(d_ij);
        }
    },
    parseRow: function(row, i){
        if (i === 0) {
            return this.parseHeader(row);
        }
        return this.parseDataRow(row);
    },
    parse: function(raw) {
        const rows = csvParseRows(raw);
        rows.forEach((r, i) => this.parseRow(r,i));
        return this.data;
    }
}

export default function dataParser(raw) {
    const dataParser = new DataParser();
    return dataParser.parse(raw);
}
