<html>
<body>
<script type="text/javascript">

$( "cite" ).each(
  function( index) {
    var name = "#Citation_" + $(this).attr("ref");
    $(this).html( "<a href=\"" + name + "\">" + $(name).text() + "</a>" );
  }
);

$( ":target" ).css( "background-color", "palegoldenrod" );

</script>
</body>
</html>
