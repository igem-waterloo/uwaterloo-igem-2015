// Implements the <cite> tag.  Requires a "ref" attribute that corresponds to one listed in references.json and that
// references.json be included in the global scope
var scrollLinks = {};

$(document).ready(function() {  
  var citation_idx = 1;
  
  // On first encounter, fill in all <cite> tags on the page with a certain reference
  $( "cite" ).each( function() {
    
    if ( $(this).html().indexOf("href") != -1 ) { // tag already been given link in reference list
    }
    else {
        var refID = $(this).attr("ref");
        
        // Add link with ordinal citation inside <cite> tag of all with refID
        $( "cite[ref~='" + refID + "']" ).each( function() {
            $(this).html( "<a class=\"ref-link\" href=\"#" + refID + "\"><em>[" + citation_idx + "]</em></a>");
        });
        
        // Add refID to list
        scrollLinks[$(this).text()] = refID;
     
        // Add appropriate li to the references list
        $("#reflist").append("<li id=\"" + refID + "\">" + references[refID] + "</li>");

        citation_idx = citation_idx + 1;
    }
  }); // end citetag.each
  $(".ref-link").each(function(i, obj) {
      $(obj).click(function() {
          scrollToAnchor(scrollLinks[obj.text]);
      });
  });
}); // end document.ready

// for smooth scrolling
function scrollToAnchor(aid){
    var aTag = $('#'+aid);
    if(aTag.length){
      $('html, body').animate({scrollTop:$(aTag).position().top}, 'slow');
    }
}
