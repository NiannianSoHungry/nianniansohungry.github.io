var posts=["2026/01/26/我目前做生信分析使用的文件树/","2026/01/27/ggplot2完全解析EP-01/"];function toRandomPost(){
    pjax.loadUrl('/'+posts[Math.floor(Math.random() * posts.length)]);
  };