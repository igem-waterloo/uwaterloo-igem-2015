$('#groupparts').load('/cgi/api/groupparts.cgi', {
    t: 'iGEM015',
    g: 'Waterloo'
}, function() {
    $('#groupparts .tablesorter').tablesorter();
});
