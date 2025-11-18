// get the ninja-keys element
const ninja = document.querySelector('ninja-keys');

// add the home and posts menu items
ninja.data = [{
    id: "nav-contact",
    title: "Contact",
    section: "Navigation",
    handler: () => {
      window.location.href = "/";
    },
  },{id: "nav-bioinformatic-analyses",
          title: "Bioinformatic analyses",
          description: "",
          section: "Navigation",
          handler: () => {
            window.location.href = "/bioinformatic-analyses/";
          },
        },{id: "nav-publications",
          title: "Publications",
          description: "",
          section: "Navigation",
          handler: () => {
            window.location.href = "/publications/";
          },
        },{id: "projects-pseudotime",
          title: 'Pseudotime',
          description: "",
          section: "Projects",handler: () => {
              window.location.href = "/projects/sc_pseudotime.html";
            },},{
      id: 'light-theme',
      title: 'Change theme to light',
      description: 'Change the theme of the site to Light',
      section: 'Theme',
      handler: () => {
        setThemeSetting("light");
      },
    },
    {
      id: 'dark-theme',
      title: 'Change theme to dark',
      description: 'Change the theme of the site to Dark',
      section: 'Theme',
      handler: () => {
        setThemeSetting("dark");
      },
    },
    {
      id: 'system-theme',
      title: 'Use system default theme',
      description: 'Change the theme of the site to System Default',
      section: 'Theme',
      handler: () => {
        setThemeSetting("system");
      },
    },];
