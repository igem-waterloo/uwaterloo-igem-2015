function showHideToggle (tag) {
    var object = $('#' + tag.id),
        icon = $('#' + tag.id + ' span');

    if ( object.hasClass('collapsed') ) {
        object.removeClass('collapsed');
        object.addClass('expanded');
        icon.removeClass('glyphicon-chevron-down');
        icon.addClass('glyphicon-chevron-up');
    } else if (object.hasClass('expanded') ) {
        object.removeClass('expanded');
        object.addClass('collapsed');
        icon.removeClass('glyphicon-chevron-up');
        icon.addClass('glyphicon-chevron-down');
    }
}
