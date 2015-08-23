// Implements the <cite> tag.  Requires a "ref" attribute that corresponds to one listed in references.json and that
// references.json be included in the global scope
$(document).ready(function() {
  $( "cite" ).each(
    function( index ) {
      var refID = $(this).attr("ref");  
    
      // Add link with ordinal citation inside <cite> tag
      $(this).html( "<a href=\"#" + refID + "\"><em>[" + (index + 1) + "]</em></a>");
    
      // Add appropriate li to the references list
      $("#reflist").append("<li id=\"" + refID + "\">" + references[refID] + "</li>")
    }
  );
});