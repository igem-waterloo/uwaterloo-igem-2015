// Implements the <cite> tag.  Requires a "ref" attribute that corresponds to one listed in references.json and that
// references.json be included in the global scope

var refLinks = {};

$(document).ready(function() {
  $( "cite" ).each(
    function( index ) {
      var refID = $(this).attr("ref");

      // Add link with ordinal citation inside <cite> tag
      $(this).html( "<a class=\"ref-link\" href=\"#" + refID + "\"><em>[" + (index + 1) + "]</em></a>");

      // Add refID to list
      scrollLinks[$(this).text()] = refID;

      // Add appropriate li to the references list
      $("#reflist").append("<li id=\"" + refID + "\">" + references[refID] + "</li>");
    }
  );

  $(".ref-link").each(function(i, obj) {
      $(obj).click(function() {
          scrollToAnchor(scrollLinks[obj.text]);
      });
  });
});

// for smooth scrolling
function scrollToAnchor(aid){
    var aTag = $('#'+aid);
    if(aTag.length){
      $('html, body').animate({scrollTop:$(aTag).position().top}, 'slow');
    }
}
