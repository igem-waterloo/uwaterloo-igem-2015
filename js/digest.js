function calculator (){
    var samples = $('.sample').map(function(i, obj){
        return obj.value;
    });
    var concentrations = $('.concentration').map(function(i, obj){
        return parseFloat(obj.value);
    });
    var amt_dna = parseFloat($('.amt_dna')[0].value),
        num_enzymes = parseInt($('.num_enzymes')[0].value),
        total_volume = parseFloat($('.total_volume')[0].value),
        table = {};

    samples.each(function(i, sample) {
        var dna = amt_dna/concentrations[i],
            enzyme = amt_dna/1000,
            buffer = total_volume/10,
            milliq = total_volume-(dna + buffer + enzyme*num_enzymes);

        table[sample] = [dna, enzyme, buffer, milliq];
    });

    samples.each(function (i, sample) {
        // check each value exists
        table[sample] = table[sample].map(function(e) {
            if (e) {
                return e.toFixed(2);
            } else {
                return 0;
            }
        });
        addColumn(sample, table[sample]);
    });
}

function addColumn (sample, col) {
    var th = document.createElement('th'),
        ids = ['res_dna', 'res_enzyme', 'res_buffer', 'res_milliq'],
        head = document.getElementById('res_label');

    head.appendChild(th);
    head.lastChild.innerHTML = sample;

    ids.forEach(function (id, i) {
        var data = document.getElementById(id),
            td = document.createElement('td');
        data.appendChild(td);
        data.lastChild.innerHTML = col[i];
    });
}

function addSample () {
    var sample = document.createElement('input'),
        concentration = document.createElement('input'),
        br = document.createElement('br'),
        samplesDiv = document.getElementById("samples");
    sample.type="text";
    sample.setAttribute("class", "sample");
    sample.setAttribute("value", "Sample Label");
    sample.style["margin-right"] = "5px";
    concentration.type="number";
    concentration.setAttribute("class", "concentration");
    concentration.setAttribute("value", "0");

    samplesDiv.appendChild(sample);
    samplesDiv.appendChild(concentration);
    samplesDiv.appendChild(br);
}
