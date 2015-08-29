// scripts for header

// detects scroll for top resizing and lower nav show
$(window).scroll(function() {
  if ($(document).scrollTop() > 50) {
    $('.main-nav').addClass('shrink');
    if ($('#inner-page-links').children().length > 0) {
      $('.navbar-lower').removeClass('hide-lower');
    }
  } else {
    $('.main-nav').removeClass('shrink');
    $('.navbar-lower').addClass('hide-lower');
  }
});

var scrollLinks = {};

// fills in lower nav with inner page links
$(document).ready(function(){
    $('.accordion-heading').addClass('link');
    $('section').addClass('link');
    $('.link').each(function(i, obj) {
        $("#inner-page-links").append('<li><a href="#" class="scroll-link">'+obj.title+'</a></li>');
        scrollLinks[obj.title] = obj.id;
    });
    $(".scroll-link").each(function(i, obj) {
        $(obj).click(function() {
            scrollToAnchor(scrollLinks[obj.text]);
        });
    });
    if ($('#inner-page-links').children().length < 1) {
        $('.navbar-lower').addClass('hide-lower');
    }
});

// for smooth scrolling
function scrollToAnchor(aid){
    var aTag = $('#'+aid);
    if(aTag.length){
      $('html, body').animate({scrollTop:$(aTag).position().top}, 'slow');
    }
}
