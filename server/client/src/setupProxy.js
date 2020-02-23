const proxy = require('http-proxy-middleware')

module.exports = function(app) {
  //app.use(proxy('/auth/google', { target: 'http://localhost:5200/' , secure: false, changeOrigin: true }));
  app.use(proxy('/api/*', { target: 'http://localhost:5200/' , secure: false, changeOrigin: true }))
  //app.use(proxy('/*.svg', { target: 'http://localhost:5000/' }))
}
