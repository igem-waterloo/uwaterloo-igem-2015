// scripts for header

// for capitalizing id
String.prototype.capitalize = function() {
    return this.charAt(0).toUpperCase() + this.slice(1);
}

// detects scroll for top resizing and lower nav show
$(window).scroll(function() {
  if ($(document).scrollTop() > 50) {
    $('.main-nav').addClass('shrink');
    if ($('#inner-page-links').children().length > 0) {
      $('.navbar-lower').removeClass('hide-lower')
    }
  } else {
    $('.main-nav').removeClass('shrink');
    $('.navbar-lower').addClass('hide-lower')
  }
});

// fills in lower nav with inner page links
$(document).ready(function(){
  $('.accordion-heading').addClass('link');
  $('.link').each(function(i, obj) {
     $("#inner-page-links").append('<li><a href="#" class="scroll-link">'+obj.id.capitalize()+'</a></li>');
  });
  $(".scroll-link").each(function(i, obj) {
    $(obj).click(function() {
      scrollToAnchor(obj.text.toLowerCase());
    });
  });
  if ($('#inner-page-links').children().length < 1) {
    $('.navbar-lower').addClass('hide-lower')
  }
});

// for smooth scrolling
function scrollToAnchor(aid){
    var aTag = $('#'+aid);
    if(aTag.length){  
      $('html, body').animate({scrollTop:$(aTag).position().top}, 'slow');
    }
}
